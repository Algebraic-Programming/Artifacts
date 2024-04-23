# Copyright 2024 Huawei Technologies Co., Ltd.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# author: Pal Andras Papp, Georg Anegg

from dataclasses import dataclass

from mip import Model, xsum, MINIMIZE, BINARY, INTEGER, OptimizationStatus

from instance import Instance
from schedule import Schedule


@dataclass
class IPSettings:
    num_supersteps: int = 3  # max number of supersteps allowed
    solver: str = 'CBC'  # CBC or GRB (if installed)
    search_emphasis: int = 0  # 2 for optimality/lower bounds, 1 for feasible solutions, 0 for balance
    cuts: int = 1  # 0 for disabling cuts (more time branching and exploring nodes), 1-3 to increase search for cuts
    mode: str = 'Full'  # 'Full' for optimizing the entire ILP, 'CS' for optimizing comm. schedule only


@dataclass
class Solver:
    instance: Instance
    ip_settings: IPSettings
    initial_solution: Schedule = None

    def __post_init__(self):
        if self.initial_solution:  # if initial solution is given, overwrite given instance with that from schedule for consistency
            self.instance = self.initial_solution.instance
        self.model = Model(solver_name=self.ip_settings.solver)
        self.model.sense = MINIMIZE
        self.model.emphasis = self.ip_settings.search_emphasis
        self.model.cuts = self.ip_settings.cuts
        self.model.round_int_vars = True
        self.model.store_search_progress_log = True
        self.model.verbose = 0
        self.total_time = 0
        self.gap = 0
        self.log = self.model.search_progress_log
        self.n = len(self.instance.G.nodes())

        self.T = self.instance.G.nodes()  # set of tasks
        self.E = self.instance.G.edges()  # set of edges
        self.P = range(self.instance.p)  # set of processors
        self.S = range(self.initial_solution.supersteps_used if self.initial_solution else self.ip_settings.num_supersteps)  # set of supersteps
        self.set_up()

    def set_up(self):
        if self.ip_settings.mode == 'CS':
            self.initializeCS()

        self.add_variables()

        if self.ip_settings.mode == 'Full':
            self.add_constraints()
        if self.ip_settings.mode == 'CS':
            self.add_CSconstraints()

        self.add_objective()
        if self.initial_solution:
            self.set_initial_solution()

    def set_initial_solution(self):
        modename = "ILP" if self.ip_settings.mode == 'Full' else "ILP(CS)"
        print(f"Starting {modename} with initial solution of (unoptimized) value {self.initial_solution.objective_value} (using {self.initial_solution.supersteps_used} supersteps)")
        if self.ip_settings.mode == 'Full':
            self.model.start = [(self.x[(task, self.initial_solution.where[task], self.initial_solution.when[task])], 1) for task in self.instance.G.nodes()]
        if self.ip_settings.mode == 'CS':
            self.model.start = [(self.cpn[(processor, self.first_needed[(task, processor)] - 1, task)], 1) for (task, processor) in self.ever_needed]

    def add_variables(self):
        # variables indicating if superstep is used at all
        self.l = {superstep: self.model.add_var(var_type=BINARY, name=f"l[{superstep}]") for superstep in self.S}

        # largest bulk communication (on any processor) at superstep
        self.c = {superstep: self.model.add_var(var_type=INTEGER, lb=0, name=f"c[{superstep}]") for superstep in self.S}

        if self.ip_settings.mode == 'Full':
            # indicator variables for tasks on processor at superstep
            self.x = {(task, processor, superstep): self.model.add_var(var_type=BINARY, name=f"x{(task, processor, superstep)}")
                      for task in self.T for processor in self.P for superstep in self.S}

            # computation load at superstep
            self.w = {superstep: self.model.add_var(var_type=INTEGER, lb=0, name=f"w[{superstep}]") for superstep in self.S}

            # communicate node i from from_processor to to_processor at superstep
            self.cppn = {(from_processor, to_processor, superstep, i): self.model.add_var(var_type=BINARY, name=f"cppn[({from_processor, to_processor, superstep, i})]")
                         for from_processor in self.P for to_processor in self.P for superstep in self.S for i in self.T}

        if self.ip_settings.mode == 'CS':
            self.cpn = {}
            self.ever_needed = set()
            for i in self.T:
                for to_processor in self.P:
                    if to_processor == self.initial_solution.where[i] or self.first_needed[(i, to_processor)] == self.initial_solution.supersteps_used + 1:
                        continue
                    self.ever_needed.add((i, to_processor))
                    for superstep in range(self.initial_solution.when[i], self.first_needed[(i, to_processor)]):
                        self.cpn[(to_processor, superstep, i)] = self.model.add_var(var_type=BINARY, name=f"cpn[({to_processor, superstep, i})]")


    def add_objective(self):
        if self.ip_settings.mode == 'Full':
            self.model.objective = xsum(self.w[superstep] + self.instance.g * self.c[superstep] + self.instance.L * self.l[superstep] for superstep in self.S) - self.instance.L
        if self.ip_settings.mode == 'CS':
            self.model.objective = xsum(self.instance.g * self.c[superstep] + self.instance.L * self.l[superstep] for superstep in self.S)

    def add_constraints(self):
        self.add_modelling_constraints()
        self.add_feasibility_constraints()
        self.add_communication_constraints()

    def add_modelling_constraints(self):
        # computational cost at superstep is max of load on processors
        for superstep in self.S:
            for processor in self.P:
                self.model.add_constr(self.w[superstep] >= xsum(self.instance.G.nodes[task]['compute_cost'] * self.x[(task, processor, superstep)] for task in self.T))

        # communication volume at superstep is max of send/receive
        for superstep in self.S:
            for processor in self.P:
                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.instance.G.nodes[task]['communication_cost'] * self.instance.numa_costs[(processor, to_processor)] * self.cppn[(processor, to_processor, superstep, task)]
                                           for task in self.T
                                           for to_processor in self.P if to_processor != processor))
                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.instance.G.nodes[task]['communication_cost'] * self.instance.numa_costs[(from_processor, processor)] * self.cppn[(from_processor, processor, superstep, task)]
                                           for task in self.T
                                           for from_processor in self.P if processor != from_processor))


    def add_feasibility_constraints(self):

        # every task is scheduled
        if self.instance.recomputation:
            for task in self.T:
                self.model.add_constr(xsum(self.x[(task, processor, superstep)] for processor in self.P for superstep in self.S) >= 1)
        else:
            for task in self.T:
                self.model.add_constr(xsum(self.x[(task, processor, superstep)] for processor in self.P for superstep in self.S) == 1)

        # use consecutive supersteps starting from 0
        for superstep in self.S[:-1]:
            self.model.add_constr(self.l[superstep] >= self.l[superstep + 1])
        self.model.add_constr(self.l[0] == 1)

        # superstep is used at all
        for superstep in self.S:
            self.model.add_constr(xsum(self.x[task, processor, superstep] for task in self.T for processor in self.P) <= self.n * self.instance.p * self.l[superstep])

        # precedence constraint: if task is computed then all of its predecessors must have been present
        for processor in self.P:
            for superstep in self.S:
                for task in self.T:
                    predecessors = list(self.instance.G.predecessors(task))
                    self.model.add_constr(
                        xsum(self.cppn[(processor, processor, superstep, predecessor)] for predecessor in predecessors) >= len(list(predecessors)) * self.x[(task, processor, superstep)])

    def add_communication_constraints(self):

        # combines two constraints: node can only be communicated if it is present; and node is present if it was computed or communicated
        for superstep in self.S:
            for processor in self.P:
                for node in self.T:
                    self.model.add_constr(self.instance.p * ((xsum(self.cppn[(from_processor, processor, superstep - 1, node)] for from_processor in self.P)
                                                              if superstep > 0 else 0) +
                                                             self.x[(node, processor, superstep)])
                                          >= xsum(self.cppn[(processor, to_processor, superstep, node)] for to_processor in self.P))

    def initializeCS(self):
        self.first_needed = {(task, processor): self.initial_solution.supersteps_used+1 for task in self.T for processor in self.P}
        for (from_node, to_node) in self.E:
            if self.initial_solution.where[from_node] == self.initial_solution.where[to_node]:
                continue
            if self.first_needed[(from_node, self.initial_solution.where[to_node])] > self.initial_solution.when[to_node]:
                self.first_needed[(from_node, self.initial_solution.where[to_node])] = self.initial_solution.when[to_node]

        self.commcost = {(task, processor): self.instance.numa_costs[(self.initial_solution.where[task], processor)] * self.instance.G.nodes[task]['communication_cost']
                         for task in self.T for processor in self.P}


    def add_CSconstraints(self):
        # communication volume at superstep is max of send/receive
        for superstep in self.S:
            for processor in self.P:

                senders = set()
                for task in self.T:
                    if (processor, superstep, task) in self.cpn.keys():
                        senders.add(task)

                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.commcost[(task, processor)] * self.cpn[(processor, superstep, task)]
                                           for task in senders))

                receivers = set()
                for task in self.T:
                    if self.initial_solution.where[task] != processor:
                        continue
                    for to_processor in self.P:
                        if (to_processor, superstep, task) in self.cpn.keys():
                            receivers.add((task, to_processor))

                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.commcost[(task, to_processor)] * self.cpn[
                                          (to_processor, superstep, task)]
                                           for (task, to_processor) in receivers))

        # every communication is scheduled
        for task in self.T:
            for processor in self.P:
                if self.first_needed[(task, processor)] == self.initial_solution.supersteps_used + 1:
                    continue
                self.model.add_constr(xsum(self.cpn[(processor, superstep, task)]
                                           for superstep in range(self.initial_solution.when[task],
                                                                  self.first_needed[(task, processor)])) == 1)

        # superstep is used at all
        for superstep in self.S:
            validsteps = set()
            for task in self.T:
                for processor in self.P:
                    if (processor, superstep, task) in self.cpn.keys():
                        validsteps.add((task, processor))

            self.model.add_constr(
                xsum(self.cpn[processor, superstep, task] for (task, processor) in validsteps) <= len(validsteps) *
                self.l[superstep])



    def optimize(self, seconds=0, minutes=0, hours=0):
        status = self.model.optimize(max_seconds=hours * 3600 + minutes * 60 + seconds)
        elapsed_time = hours * 3600 + minutes * 60 + seconds
        gap = 99.9999

        if status == OptimizationStatus.OPTIMAL:
            gap = 0
            elapsed_time = self.log.log[-1][0] if self.log.log else 0

        if status == OptimizationStatus.FEASIBLE:
            if self.log.log:
                lb, ub = self.log.log[-1][1]
                gap = abs(ub - lb) / ub
        if status == OptimizationStatus.NO_SOLUTION_FOUND:
            print("Not enough time to set up. Increase time given for optimization step.")
        self.total_time += elapsed_time
        self.gap = gap
        addednote = ' for CS' if self.ip_settings.mode == "CS" else ''
        supsteps = max(step for step in self.S if self.l[step].x == 1) + (2 if self.ip_settings.mode == "CS" else 1)
        print(f"    Best solution: objective value {round(self.model.objective_value) if self.model.objective_value else -1:>5} ({status._name_:>8}{addednote}, gap {100 * self.model.gap:>5.2f}%) after {elapsed_time:>6.2f} seconds  ({supsteps} supersteps)")
        return status, elapsed_time, gap

    def get_schedule(self):
        copied_dag = self.instance.G.copy()
        copied_instance = Instance(copied_dag, self.instance.L, self.instance.g, self.instance.p, self.instance.numa_costs)

        if self.ip_settings.mode == 'Full':
            proc_data = {node: next(processor for processor in self.P for superstep in self.S
                                if self.x[(node, processor, superstep)].x == 1) for node in self.instance.G.nodes()}
            superstep_data = {node: next(superstep for processor in self.P for superstep in self.S
                                if self.x[(node, processor, superstep)].x == 1) for node in self.instance.G.nodes()}
        if self.ip_settings.mode == 'CS':
            proc_data = {node: self.initial_solution.where[node] for node in self.instance.G.nodes()}
            superstep_data = {node: self.initial_solution.when[node] for node in self.instance.G.nodes()}

        comm_schedule = {}

        if self.ip_settings.mode == 'Full':
            for from_proc in self.P:
                for to_proc in self.P:
                    if from_proc == to_proc:
                        continue
                    for superstep in self.S:
                        for task in self.T:
                            if self.cppn[(from_proc, to_proc, superstep, task)].x == 1:
                                comm_schedule[(task, from_proc, to_proc)]=superstep
        if self.ip_settings.mode == 'CS':
            for task in self.T:
                for superstep in self.S:
                    for to_proc in self.P:
                        if to_proc == self.initial_solution.where[task] or (
                        to_proc, superstep, task) not in self.cpn.keys():
                            continue
                        if self.cpn[(to_proc, superstep, task)].x == 1:
                            comm_schedule[(task, self.initial_solution.where[task], to_proc)] = superstep

        return Schedule(copied_instance, proc_data, superstep_data, comm_schedule)


def initialized_optimization(initial_solution,
                             ip_settings,
                             total_seconds=0, total_minutes=0, total_hours=0,
                             ):
    total_time = 3600 * total_hours + 60 * total_minutes + total_seconds
    ip_settings = IPSettings(solver=ip_settings.solver,
                             num_supersteps=initial_solution.supersteps_used,
                             search_emphasis=ip_settings.search_emphasis,
                             mode=ip_settings.mode)
    solver = Solver(initial_solution.instance, ip_settings, initial_solution)
    print(f'Number of variables: {solver.model.num_cols}')
    solver.optimize(seconds=total_time)
    return solver


