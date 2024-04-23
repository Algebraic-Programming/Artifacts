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
from typing import Dict
from shutil import copy
from os import path

import networkx as nx

import instance
from instance import Instance, read_instance


@dataclass
class Schedule:
    instance: Instance  # the DAG in instance.G has node attributes processor and superstep
    where: Dict[int, int]
    when: Dict[int, int]
    communication_schedule: Dict = None  # the communication schedule with entries (node_id, from_processor, to_processor): superstep

    def __post_init__(self):
        self.supersteps_used = max(self.when.values()) + 1
        self.T = self.instance.G.nodes()
        self.P = range(self.instance.p)
        self.E = self.instance.G.edges()
        self.S = range(self.supersteps_used)

        self.invalid = False
        for node in self.T:
            if self.where[node] == -1:
                self.invalid = True

        if not self.invalid:
            if self.communication_schedule == None:
                self.communication_schedule = self.lazy_communication_schedule()
            self.objective_value = self.compute_objective()


    # checks the compute schedule for validity
    def validate(self):
        for source, target in self.instance.G.edges():
            required_diff = 0 if self.where[source] == self.where[target] else 1
            if self.when[target] < self.when[source] + required_diff:
                print(f'Error in edge {source}-{target}: supersteps are {self.when[source]} {self.when[target]}, processors are {self.where[source]} {self.where[target]}')
                return False

        if self.communication_schedule is not None:
            first_at = {(task, processor): self.supersteps_used for task in self.T for processor in self.P}
            for task in self.T:
                first_at[(task, self.where[task])] = self.when[task]
            for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
                first_at[(task, to_proc)] = min(first_at[(task, to_proc)], superstep + 1)

            for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
                if superstep < first_at[(task, from_proc)]:
                    print(f'Error in edge {(task, from_proc, to_proc)}-{superstep} with {first_at[(task, from_proc)]}')
                    return False
            for source, target in self.instance.G.edges():
                if first_at[(source, self.where[target])] > self.when[target]:
                    print(f'Error in edge {(source, target)}: {source} is not yet present on proc {self.where[target]} in step {self.when[target]}')
                    return False

        return True

    def lazy_communication_schedule(self):
        self.temp_schedule = {}
        for from_node, to_node in self.E:
            if self.where[from_node] == self.where[to_node]:
                continue
            if self.temp_schedule.get((from_node, self.where[from_node], self.where[to_node])) is not None:
                self.temp_schedule[(from_node, self.where[from_node], self.where[to_node])] = min(self.when[to_node]-1, self.temp_schedule[(from_node, self.where[from_node], self.where[to_node])])
            else:
                self.temp_schedule[(from_node, self.where[from_node], self.where[to_node])] = self.when[to_node]-1
        return self.temp_schedule

    def set_lazy_schedule(self):
        self.communication_schedule = self.lazy_communication_schedule()

    def compute_objective(self):
        if self.invalid:
            return -1
        self.cpps = {(from_proc, to_proc, superstep): 0 for from_proc in self.P for to_proc in self.P for superstep in self.S}
        for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
            self.cpps[(from_proc, to_proc, superstep)]+=self.instance.numa_costs[from_proc, to_proc]*self.instance.G.nodes()[task]['communication_cost']
        self.rec = {(processor, superstep): sum(self.cpps[(from_processor, processor, superstep)] for from_processor in self.P if from_processor != processor) for processor in self.P for superstep in
                    self.S}
        self.send = {(processor, superstep): sum(self.cpps[(processor, to_processor, superstep)] for to_processor in self.P if to_processor != processor) for processor in self.P for superstep in
                     self.S}
        self.C = {superstep: max(max(self.rec[(processor, superstep)], self.send[(processor, superstep)]) for processor in self.P) for superstep in self.S}
        self.L = {superstep: 1 if self.C[superstep]>0 else 0 for superstep in self.S}
        self.wp = {(processor, superstep): 0 for processor in self.P for superstep in self.S}
        for task in self.T:
            self.wp[(self.where[task], self.when[task])] += self.instance.G.nodes()[task]['compute_cost']
        self.W = {superstep: max(self.wp[(processor, superstep)] for processor in self.P) for superstep in self.S}
        # uncomment for details on cost
        #Wsum = sum(self.W[superstep] for superstep in self.S)
        #Csum = sum(self.instance.g * self.C[superstep] for superstep in self.S)
        #Lsum = sum(self.instance.L * self.L[superstep] for superstep in self.S)
        #print(f'Work: {Wsum}')
        #print(f'Comm: {Csum}')
        #print(f'Latency: {Lsum}')
        #print(f'Total: {Wsum + Csum + Lsum}')
        return sum(self.W[superstep] + self.instance.g * self.C[superstep] + self.instance.L * self.L[superstep] for superstep in self.S)

    def remove_needless_supersteps(self):
        self.L = {superstep: False for superstep in self.S}
        for superstep in self.communication_schedule.values():
            self.L[superstep] = True
        newIndex = {superstep: 0 for superstep in self.S}
        currentIdx = 0
        for superstep in self.S:
            newIndex[superstep] = currentIdx
            if self.L[superstep]:
                currentIdx += 1

        for node in self.T:
            self.when[node] = newIndex[self.when[node]]
        for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
            self.communication_schedule[(task, from_proc, to_proc)] = newIndex[superstep]
        self.supersteps_used = max(self.when.values()) + 1
        self.S = range(self.supersteps_used)
        self.objective_value = self.compute_objective()

    def lazy_fill_comm_schedule(self):
        if self.communication_schedule is None or self.invalid:
            return

        first_at = {(task, processor): self.supersteps_used for task in self.T for processor in self.P}
        for task in self.T:
            first_at[(task, self.where[task])] = self.when[task]
        for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
            first_at[(task, to_proc)] = min(first_at[(task, to_proc)], superstep + 1)

        first_needed = {(task, processor): self.supersteps_used for task in self.T for processor in self.P}
        for source, target in self.instance.G.edges():
            first_needed[(source, self.where[target])] = min(first_needed[(source, self.where[target])], self.when[target])

        for source, target in self.instance.G.edges():
            if first_at[(source, self.where[target])] > first_needed[(source, self.where[target])]:
                self.communication_schedule[(source, self.where[source], self.where[target])] = first_needed[(source, self.where[target])] - 1

    def clean_comm_schedule(self):
        if self.communication_schedule is None or self.invalid:
            return

        # remove duplicate sending
        for node in self.T:
            for from_proc in self.P:
                if self.communication_schedule.get((node, from_proc, self.where[node])) is not None:
                    self.communication_schedule.pop((node, from_proc, self.where[node]))

        for (node, from_proc, to_proc), superstep in self.communication_schedule.items():
            if from_proc == to_proc:
                print("Warning: unneeded communication step.")

        in_edges = {node: [] for node in self.T}
        for source, target in self.instance.G.edges():
            in_edges[target].append(source)
        comm_by_node_proc = {(node, proc): set() for node in self.T for proc in self.P}
        for (node, from_proc, to_proc), superstep in self.communication_schedule.items():
            comm_by_node_proc[(node, to_proc)].add((from_proc, superstep))
        for node in self.T:
            for to_proc in self.P:
                if len(comm_by_node_proc[(node, to_proc)]) > 1:
                    best_from = -1
                    best_step = -1
                    for from_proc, superstep in comm_by_node_proc[(node, to_proc)]:
                        if best_step == -1 or superstep < best_step or (superstep == best_step and self.instance.numa_costs[(from_proc, to_proc)] < self.instance.numa_costs[(best_from, to_proc)]):
                            best_from = from_proc
                            best_step = superstep
                    for from_proc, superstep in comm_by_node_proc[(node, to_proc)]:
                        if from_proc != best_from or superstep != best_step:
                            self.communication_schedule.pop((node, from_proc, to_proc))

        # remove unnecessary comm steps
        comm_by_superstep = {superstep: set() for superstep in self.S}
        for (node, from_proc, to_proc), superstep in self.communication_schedule.items():
            comm_by_superstep[superstep].add((node, from_proc, to_proc))
        node_by_superstep = {superstep: set() for superstep in self.S}
        for node in self.T:
            node_by_superstep[self.when[node]].add(node)

        needed = {(node, proc): False for node in self.T for proc in self.P}
        for superstep in range(self.supersteps_used-1, -1, -1):
            for (node, from_proc, to_proc) in comm_by_superstep[superstep]:
                if not needed[(node, to_proc)]:
                    self.communication_schedule.pop((node, from_proc, to_proc))
            for node in node_by_superstep[superstep]:
                for pred in in_edges[node]:
                    needed[(pred, self.where[node])] = True
            for (node, from_proc, _) in comm_by_superstep[superstep]:
                needed[(node, from_proc)] = True


    def write_file(self, input_filename, output_filename):
        copy(input_filename, output_filename)
        node_data = [f"{task} {processor} {superstep}"
                     for task in self.T for processor in self.P for superstep in self.S
                     if self.where[task] == processor
                     and self.when[task] == superstep]
        communication_data = []
        if self.communication_schedule is not None:
            communication_data.append(f"{len(self.communication_schedule)}")
            for (task, from_proc, to_proc), superstep in self.communication_schedule.items():
                communication_data.append(f"{task} {from_proc} {to_proc} {superstep}")
        with open(output_filename, 'a') as file:
            file.write('% BSP Schedule Data\n')
            file.write('\n'.join(node_data))
            file.write('\n')
            if self.communication_schedule is not None:
                file.write('% Communication Schedule\n')
                file.write('\n'.join(communication_data))
                file.write('\n')
        print(f"Output written to {output_filename}")

    def CanOptimizeComm(self):
        return self.supersteps_used >= 3

    def GetNodesInSupersteps(self, fromSupstep, toSupstep):
        nodes = []
        for node in self.instance.G.nodes():
            if fromSupstep <= self.when[node] <= toSupstep:
                nodes.append(node)
        return nodes

# SolutionList: stores the schedules found for a given instance
@dataclass
class SolutionList:
    instance: Instance  # problem
    schedules: Dict[str, Schedule]

    def __post_init__(self):
        self.lower_bound = 0
        if self.instance is not None:
            self.AddTrivial()

    def add(self, name, sched):
        self.schedules[name] = sched

    def AddTrivial(self):
        processor = {v: 0 for v in self.instance.G.nodes()}
        superstep = {v: 0 for v in self.instance.G.nodes()}
        self.schedules["Trivial"] = Schedule(self.instance, processor, superstep, None)

    def GetBestSchedule(self):
        best = None
        for algo, sched in self.schedules.items():
            if algo != "Trivial" and algo != "cilk" and algo != "ETF" and algo != "BL-EST" and (best is None or sched.objective_value < best.objective_value):
                best = sched
        if best is None:
            best = self.schedules["Trivial"]
        return best

    def GetBestName(self):
        best = None
        name = None
        for algo, sched in self.schedules.items():
            if algo != "Trivial" and algo != "cilk" and algo != "ETF" and algo != "BL-EST" and (best is None or sched.objective_value < best.objective_value):
                best = sched
                name = algo
        return name

    def Print(self):
        for (algo, sched) in self.schedules.items():
            if sched.validate() is True:
                print(f'{algo}: {sched.objective_value} (in {sched.supersteps_used} supersteps) - {len(sched.communication_schedule)/(sched.supersteps_used-1 if sched.supersteps_used > 1 else 1)}')
            else:
                print(f'ERROR: schedule returned by {algo} is invalid.')
        print(f'Best lower bound: {self.lower_bound}')

    def PrintCSV(self, instance_name, filename):
        file = open(filename, 'a+')
        file.write(f'\n{instance_name}\n')
        for (algo, sched) in self.schedules.items():
            file.write(f'{algo},{sched.objective_value},{sched.supersteps_used}\n')
        file.write(f'{self.lower_bound}\n')
        file.close()

# SolutionDict: dictionary of instance names to SolutionLists
def CreateSolutionDict(instance_names, settings):
    sol = {}
    for instance_path in instance_names:
        name = path.split(instance_path)[1]
        sol[name] = SolutionList(read_instance(instance_path, settings.NUMA), {})
    return sol


def read_schedule(filename, NUMA=True, validate=False):
    cleaned = instance.GetCleanedFile(filename)

    inst = instance.GetInstanceFromCleanedFile(cleaned, NUMA)

    num_hyperedges, num_nodes, num_pins = cleaned[0]
    p = cleaned[num_pins + num_nodes + 1][0]
    numa_lines = p * p if NUMA else 0
    schedule_start = num_pins + num_nodes + 2 + numa_lines

    schedule_info = cleaned[schedule_start:schedule_start + num_nodes]
    processors = {entry[0]: entry[1] for entry in schedule_info}
    supersteps = {entry[0]: entry[2] for entry in schedule_info}

    try:  # check if schedule does not have comm schedule info
        num_communications = cleaned[schedule_start + num_nodes][0]
        comm_schedule_info = cleaned[schedule_start + num_nodes + 1:schedule_start + num_nodes + num_communications+1]
        communication_schedule = {(node, from_processor, to_processor): superstep for node, from_processor, to_processor, superstep in comm_schedule_info}
    except:
        communication_schedule = None

    s = Schedule(inst, processors, supersteps, communication_schedule)
    if validate:
        if not s.validate():
            print(f"Schedule {filename} is not valid.")
    return s

# auxiliary for multilevel
def get_uncontracted_schedule(inst, schedule, contraction_file):
    cleaned = instance.GetCleanedFile(contraction_file)
    old_n, new_n = cleaned[0]
    image_lines = cleaned[old_n-new_n+1:2*old_n-new_n+1]
    image = {}
    for line in image_lines:
        image[line[0]] = line[1]

    proc_data = {}
    superstep_data = {}
    for node in inst.G.nodes():
        proc_data[node] = schedule.where[image[node]]
        superstep_data[node] = schedule.when[image[node]]

    if schedule.communication_schedule is None:
        return Schedule(inst, proc_data, superstep_data, None)

    contains = {node: set() for node in schedule.T}
    for node in inst.G.nodes():
        contains[image[node]].add(node)
    comm_schedule = {}
    for (target, from_proc, to_proc), superstep in schedule.communication_schedule.items():
        for node in contains[target]:
            comm_schedule[(node, from_proc, to_proc)] = superstep
    new_schedule = Schedule(inst, proc_data, superstep_data, comm_schedule)
    new_schedule.clean_comm_schedule()
    return new_schedule


def get_uncontracted_cost(instance, schedule, contraction_file):
    schedule = get_uncontracted_schedule(instance, schedule, contraction_file)
    return schedule.compute_objective()
