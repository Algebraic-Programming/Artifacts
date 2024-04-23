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

import copy
from schedule import Schedule


@dataclass
class PartIPSettings:
    num_supersteps: int = 3  # max number of supersteps allowed
    solver: str = 'CBC'  # CBC or GRB (if installed)
    search_emphasis: int = 2  # 2 for optimality/lower bounds, 1 for feasible solutions, 0 for balance
    cuts: int = 1  # 0 for disabling cuts (more time branching and exploring nodes), 1-3 to increase search for cuts


@dataclass
class PartSolver:
    ip_settings: PartIPSettings
    nodes: []  # subset of nodes to schedule in the ILP

    fromSupstep: int = 0
    toSupstep: int = 1
    schedule: Schedule = None

    def __post_init__(self):
        self.free = {node: False for node in self.schedule.T}
        for node in self.nodes:
            self.free[node] = True

        self.S = range(self.fromSupstep, self.toSupstep + 1)
        self.S_plus = range(max(0, self.fromSupstep-1), self.toSupstep + 1) # S, extended with the communication phase right before S

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
        self.status = None

        self.instance = self.schedule.instance

        self.n = len(self.instance.G.nodes())
        self.T = self.schedule.T  # set of tasks
        self.E = self.schedule.instance.G.edges()  # set of edges
        self.P = self.schedule.P  # set of processors

        self.compute_auxiliary_graph_data()
        self.schedule.compute_objective()

        self.direct_only_simple = False # force nodes to be sent always from where they were computed - currently unused
        self.set_up()

    def compute_auxiliary_graph_data(self):
        self.out_edges = {node: [] for node in self.instance.G.nodes()}
        self.in_edges = {node: [] for node in self.instance.G.nodes()}
        for source, target in self.schedule.instance.G.edges():
            self.out_edges[source].append(target)
            self.in_edges[target].append(source)

        self.predecessors = set()
        for node in self.nodes:
            for pred in self.in_edges[node]:
                if not self.free[pred] and self.schedule.where[pred] != -1:
                    self.predecessors.add(pred)

    def set_up(self):

        self.add_variables()
        self.add_constraints()
        self.add_objective()

        self.all_set = True
        for node in self.nodes:
            if self.schedule.where[node] == -1:
                self.all_set = False
        if self.all_set:
            self.set_initial_solution()


    def set_initial_solution(self):
        print(f"Starting ILP with initial solution of (unoptimized) value {self.schedule.objective_value} (using {self.schedule.supersteps_used} supersteps)")
        self.model.start = [(self.x[(node, self.schedule.where[node], self.schedule.when[node])], 1) for node in self.nodes]

    def add_variables(self):
        # variables indicating if superstep is used at all
        self.l = {superstep: self.model.add_var(var_type=BINARY, name=f"l[{superstep}]") for superstep in self.S_plus}

        # largest bulk communication (on any processor) at superstep
        self.c = {superstep: self.model.add_var(var_type=INTEGER, lb=0, name=f"c[{superstep}]") for superstep in self.S_plus}

        # indicator variables for tasks on processor at superstep
        self.x = {(task, processor, superstep): self.model.add_var(var_type=BINARY, name=f"x{(task, processor, superstep)}")
                  for task in self.nodes for processor in self.P for superstep in self.S}

        # computation load at superstep
        self.w = {superstep: self.model.add_var(var_type=INTEGER, lb=0, name=f"w[{superstep}]") for superstep in self.S}

        # communication options for the (movable) nodes
        self.cppn = {}
        self.needed_after = set()
        for node in self.nodes:
            has_free_successor = False
            successors_after = {proc: False for proc in self.P}
            for succ in self.out_edges[node]:
                if self.free[succ]:
                    has_free_successor = True
                elif self.schedule.where[succ] != -1:
                    successors_after[self.schedule.where[succ]] = True

            if has_free_successor:
                for from_proc in self.P:
                    for to_proc in self.P:
                        for superstep in range(self.fromSupstep, self.toSupstep):
                            self.cppn[(from_proc, to_proc, superstep, node)] = self.model.add_var(var_type=BINARY, name=f"cppn[({from_proc, to_proc, superstep, node})]")

            for proc in self.P:
                self.cppn[(proc, proc, self.toSupstep, node)] = self.model.add_var(var_type=BINARY,name=f"cppn[({proc, proc, self.toSupstep, node})]")

            for proc in self.P:
                if not successors_after[proc]:
                    continue
                self.needed_after.add((node, proc))
                for from_proc in self.P:
                    if self.cppn.get((from_proc, proc, self.toSupstep, node)) is None:
                        self.cppn[(from_proc, proc, self.toSupstep, node)] = self.model.add_var(var_type=BINARY, name=f"cppn[({from_proc, proc, self.toSupstep, node})]")

        # communication options for predecessors of the (movable) nodes
        present_before = {(pred, proc): False for pred in self.predecessors for proc in self.P}
        successor_on_proc = {(pred, proc): False for pred in self.predecessors for proc in self.P}
        originally_sent = {(pred, proc): False for pred in self.predecessors for proc in self.P}
        for pred in self.predecessors:
            present_before[(pred, self.schedule.where[pred])] = True
            for succ in self.out_edges[pred]:
                if self.schedule.when[succ] > self.toSupstep:
                    successor_on_proc[(pred, self.schedule.where[succ])] = True
        for (node, from_proc, to_proc), superstep in self.schedule.communication_schedule.items():
            if node in self.predecessors:
                if superstep < self.fromSupstep:
                    present_before[(node, to_proc)] = True
                elif superstep <= self.toSupstep:
                    originally_sent[(node, to_proc)] = True
        for pred in self.predecessors:
            for proc in self.P:
                if present_before[(pred, proc)]:
                    continue
                for superstep in self.S_plus:
                    self.cppn[(self.schedule.where[pred], proc, superstep, pred)] = self.model.add_var(var_type=BINARY, name=f"cppn[({self.schedule.where[pred], proc, superstep, pred})]")
                for superstep in self.S:
                    self.cppn[(proc, proc, superstep, pred)] = self.model.add_var(var_type=BINARY, name=f"cppn[({proc, proc, superstep, pred})]")
                if successor_on_proc[(pred, proc)] and originally_sent[(pred, proc)]:
                    self.needed_after.add((pred, proc))


    def add_objective(self):
        self.model.objective = xsum(self.w[superstep] for superstep in self.S) + xsum(self.instance.g * self.c[superstep] + self.instance.L * self.l[superstep] for superstep in self.S_plus)

    def add_constraints(self):
        self.add_modelling_constraints()
        self.add_feasibility_constraints()
        self.add_communication_constraints()

    def add_modelling_constraints(self):
        # computational cost at superstep is max of load on processors
        for superstep in self.S:
            for processor in self.P:
                self.model.add_constr(self.w[superstep] >= xsum(self.instance.G.nodes[node]['compute_cost'] * self.x[(node, processor, superstep)] for node in self.nodes))

        # auxiliary for fixed comm costs
        fixed_send_cost = {(proc, superstep): 0 for proc in self.P for superstep in self.S_plus}
        fixed_rec_cost = {(proc, superstep): 0 for proc in self.P for superstep in self.S_plus}
        for (node, from_proc, to_proc), superstep in self.schedule.communication_schedule.items():
            if self.free[node] or node in self.predecessors or superstep not in self.S_plus:
                continue
            fixed_send_cost[(from_proc, superstep)] += self.instance.G.nodes[node]['communication_cost'] * self.instance.numa_costs[(from_proc, to_proc)]
            fixed_rec_cost[(to_proc, superstep)] += self.instance.G.nodes[node]['communication_cost'] * self.instance.numa_costs[(from_proc, to_proc)]

        # storing cppn variables
        cppn_sent = {(proc, superstep): [] for proc in self.P for superstep in self.S_plus}
        cppn_rec = {(proc, superstep): [] for proc in self.P for superstep in self.S_plus}
        for (from_proc, to_proc, superstep, node) in self.cppn.keys():
            if from_proc != to_proc:
                cppn_sent[(from_proc, superstep)].append((to_proc, node))
                cppn_rec[(to_proc, superstep)].append((from_proc, node))


        # communication volume at superstep is max of send/receive
        for superstep in self.S_plus:
            for processor in self.P:
                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.instance.G.nodes[node]['communication_cost'] * self.instance.numa_costs[(processor, to_processor)] * self.cppn[(processor, to_processor, superstep, node)]
                                           for (to_processor, node) in cppn_sent[(processor, superstep)])
                                      + fixed_send_cost[(processor, superstep)])

                self.model.add_constr(self.c[superstep] >=
                                      xsum(self.instance.G.nodes[node]['communication_cost'] * self.instance.numa_costs[(from_processor, processor)] * self.cppn[(from_processor, processor, superstep, node)]
                                           for (from_processor, node) in cppn_rec[(processor, superstep)])
                                      + fixed_rec_cost[(processor, superstep)])

        # for direct only mode (currently unused)
        if self.direct_only_simple:
            for (from_proc, to_proc, superstep, node) in self.cppn.keys():
                if from_proc == to_proc or not self.free[node]:
                    continue
                self.model.add_constr(self.cppn[(from_proc, to_proc, superstep, node)] <= xsum(self.x[(node, from_proc, superstep)] for superstep in self.S))


    def add_feasibility_constraints(self):

        # every task is scheduled
        for node in self.nodes:
            self.model.add_constr(xsum(self.x[(node, processor, superstep)] for processor in self.P for superstep in self.S) == 1)

        #superstep is used at all
        max_numa = max(self.instance.numa_costs[(from_proc, to_proc)] for from_proc in self.P for to_proc in self.P)
        max_comm_weight = max(self.instance.G.nodes[node]['communication_cost'] for node in self.instance.G.nodes())
        for superstep in self.S_plus:
            self.model.add_constr(self.c[superstep] <= len(self.instance.G.nodes()) * max_numa * max_comm_weight * self.instance.p * self.l[superstep])

        # precedence constraint: if task is computed then all of its predecessors must have been present
        for node in self.nodes:
            for pred in self.in_edges[node]:
                for processor in self.P:
                    for superstep in self.S:
                        if self.cppn.get((processor, processor, superstep, pred)) is not None:
                            self.model.add_constr(self.cppn[(processor, processor, superstep, pred)] >= self.x[(node, processor, superstep)])


    def add_communication_constraints(self):

        # general constraints -  combines two constraints: node can only be communicated if it is present; and node is present if it was computed or communicated
        # [the fact that we need to add these constraints is indicated by the fact that (proc, poc, superstesp, node) was added to cppn before]
        for (temp_proc, proc, superstep, node) in self.cppn.keys():
            if temp_proc != proc:
                continue
            
            source_options = []
            for from_proc in self.P:
                if self.cppn.get((from_proc, proc, superstep-1, node)) is not None:
                    source_options.append(from_proc)

            for to_proc in self.P:
                if self.cppn.get((proc, to_proc, superstep, node)) is not None:
                    self.model.add_constr((xsum(self.cppn[(from_proc, proc, superstep - 1, node)] for from_proc in source_options) if len(source_options) > 0 else 0) +
                        (self.x[(node, proc, superstep)] if self.free[node] else 0)
                        >= self.cppn[(proc, to_proc, superstep, node)])

        # constraints on last superstep
        # first: check is the values needed afterwards are actually sent later (after toSupstep)
        sent_after = {(node, from_proc, to_proc): False for node in self.nodes for from_proc in self.P for to_proc in self.P}
        for (node, from_proc, to_proc), superstep in self.schedule.communication_schedule.items():
            if self.free[node] and superstep > self.toSupstep:
                sent_after[(node, from_proc, to_proc)] = True

        # add constraints on last superstep
        for node, proc in self.needed_after:
                source_options = []
                for from_proc in self.P:
                    if self.cppn.get((from_proc, proc, self.toSupstep, node)) is not None:
                        source_options.append(from_proc)
                if self.free[node]:
                    new_options = [from_proc for from_proc in self.P if not sent_after[(node, from_proc, proc)]]
                    self.model.add_constr(
                        xsum(self.cppn[(from_proc, proc, self.toSupstep, node)] for from_proc in source_options) >= xsum(self.x[(node, from_proc, superstep)] for from_proc in new_options for superstep in self.S))
                else:
                    self.model.add_constr(xsum(self.cppn[(from_proc, proc, self.toSupstep, node)] for from_proc in source_options) >= 1)



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
        self.status = status
        self.total_time += elapsed_time
        self.gap = gap
        supsteps = 0
        for superstep in self.S:
            if self.l[superstep].x == 1:
                supsteps += 1
        print(f"    Best solution: objective value {round(self.model.objective_value) if self.model.objective_value else -1:>5} ({status._name_:>8}, gap {100 * self.model.gap:>5.2f}%) after {elapsed_time:>6.2f} seconds  ({supsteps} supersteps)")
        return status, elapsed_time, gap

    def get_schedule(self):

        if self.status == OptimizationStatus.NO_SOLUTION_FOUND:
            return None

        proc_data = {node: self.schedule.where[node] for node in self.instance.G.nodes()}
        superstep_data = {node: self.schedule.when[node] for node in self.instance.G.nodes()}

        for node in self.instance.G.nodes():
            if self.free[node]:
                for proc in self.P:
                    for superstep in self.S:
                        if self.x[(node, proc, superstep)].x == 1:
                            proc_data[node] = proc
                            superstep_data[node] = superstep

        comm_schedule = copy.deepcopy(self.schedule.communication_schedule)
        # cleaning comm. schedule
        if self.all_set:
            to_delete = []
            for (node, from_proc, to_proc), superstep in comm_schedule.items():
                # reassign free nodes completely
                if self.free[node]:
                    to_delete.append((node, from_proc, to_proc))
                # more intricate
                if node not in self.predecessors:
                    continue
                unneeded = True
                for succ in self.out_edges[node]:
                    if proc_data[succ] == to_proc:
                        unneeded = False
                for other_node, other_from_proc, _ in comm_schedule.keys():
                    if other_node == node and other_from_proc == to_proc:
                        unneeded = False
                if unneeded:
                    to_delete.append((node, from_proc, to_proc))
            for entry in to_delete:
                comm_schedule.pop(entry)

        # the reason why this part cannot be done in a simpler way: the ILP solver sometimes
        # communicates the same value twice between processors when this is not suboptimal
        comm_steps = {}
        for (from_proc, to_proc, superstep, node) in self.cppn.keys():
            if from_proc != to_proc and self.cppn[(from_proc, to_proc, superstep, node)].x == 1:
                if comm_steps.get((node, from_proc, to_proc)) is None:
                    comm_steps[(node, from_proc, to_proc)] = [superstep]
                else:
                    comm_steps[(node, from_proc, to_proc)].append(superstep)

        for (node, from_proc, to_proc) in comm_steps.keys():
            comm_schedule[(node, from_proc, to_proc)] = min(comm_steps[(node, from_proc, to_proc)])

        new_schedule = Schedule(self.instance, proc_data, superstep_data, comm_schedule)
        if not new_schedule.invalid:
            new_schedule.lazy_fill_comm_schedule()
            new_schedule.clean_comm_schedule()
            if not new_schedule.validate():
                print("Error filling new comm schedule.")
        new_schedule.objective_value = new_schedule.compute_objective()
        return new_schedule


def initialized_optimization_partial(initial_solution,
                             ip_settings,
                             nodeSet, fromSupstep, toSupstep,
                             total_seconds=0, total_minutes=0, total_hours=0,
                             ):
    total_time = 3600 * total_hours + 60 * total_minutes + total_seconds
    ip_settings = PartIPSettings(solver=ip_settings.solver,
                             num_supersteps=initial_solution.supersteps_used,
                             search_emphasis=ip_settings.search_emphasis)
    solver = PartSolver(ip_settings, nodeSet, fromSupstep, toSupstep, initial_solution)
    print(f'Number of variables: {solver.model.num_cols}')
    solver.optimize(seconds=total_time)
    return solver


