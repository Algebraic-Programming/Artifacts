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

from milppart import PartIPSettings, initialized_optimization_partial
from milp import IPSettings, initialized_optimization
from schedule import read_schedule, SolutionList, Schedule
from instance import read_instance
from dataclasses import dataclass
import multiprocessing
import queue
import util
import os
import copy
from math import ceil
from sys import platform

per = "/" if platform == "linux" else "\\"


# parameters for partial ILP representation
@dataclass
class PartialParams:
    nodeSet: []
    fromStep: int = 0
    toStep: int = 0


# main functions for solving an ILP problem
def solving_ILP(mode, initial, time_limit_min, partial_params, queue):
    if mode == "Partial":
        ip_settings = PartIPSettings(search_emphasis=2, cuts=3)
        solver = initialized_optimization_partial(
            initial_solution=initial,
            ip_settings=ip_settings,
            nodeSet=partial_params.nodeSet,
            fromSupstep=partial_params.fromStep,
            toSupstep=partial_params.toStep,
            total_minutes=time_limit_min
        )
    else: #full or CS
        ip_settings = IPSettings(search_emphasis=2, cuts=3, mode=mode)
        solver = initialized_optimization(
            initial_solution=initial,
            ip_settings=ip_settings,
            total_minutes=time_limit_min
        )

    queue.put((solver.get_schedule(), solver.model.objective_bound))


def time_limited_ILP(mode, initial_sol, time_limit_min, partial_params):
    Q = multiprocessing.Queue()
    process = multiprocessing.Process(target=solving_ILP, args=(mode, initial_sol, time_limit_min, partial_params, Q))
    process.start()
    try:
        output = Q.get(timeout=time_limit_min * 60 + 10)
    except queue.Empty:
        output = (None, 0)

    if process.is_alive():
        process.terminate()
        process.join()

    if output[0] is not None:
        return (True,) + output
    else:
        return (False,) + output


# wrapper for solving a single ILP
def run_ILP(name, mode, initial_sol, time_limit_min=1, partial_params = None):

    try:
        success, s, lower_bound = time_limited_ILP(mode, initial_sol, time_limit_min, partial_params)
        if success:
            return s, lower_bound
        else:
            print(f'Time out during ILP solving of instance {name}. Force-shutdown solver.')
            return None, 0

    except BaseException as e:
        print(f"Error in {name}: {e}")


# ILPfull
def ILP_full_single(instance_name, solDict, timeLimit):
    name = os.path.split(instance_name)[1]
    # find best initial solution
    best_initial = solDict[name].GetBestSchedule()

    # run full ILP
    (schedule, lower_bound) = run_ILP(name, "Full", best_initial, timeLimit)
    if schedule is not None:
        solDict[name].add('ILP', schedule)
        solDict[name].lower_bound = lower_bound
        schedule.write_file(f'{instance_name}', f'..{per}schedules{per}ILP_{name}')


# ILPcs
def ILP_CS_single(instance_name, solDict, timeLimit):
    name = os.path.split(instance_name)[1]
    # find best initial solution
    best_initial = solDict[name].GetBestSchedule()
    if best_initial.CanOptimizeComm():
        # run ILP for CS
        schedule, _ = run_ILP(name, "CS", best_initial, timeLimit)
        if schedule is not None:
            solDict[name].add('ILPCS', schedule)
            schedule.write_file(f'{instance_name}', f'..{per}schedules{per}ILPCS_{name}')


# ILPpart
def ILP_iter_single(instance_name, solDict, timeLimit, batch_size=40):
    name = os.path.split(instance_name)[1]
    # find best initial solution
    schedule = solDict[name].GetBestSchedule()

    nodesInSupstep = {supstep: 0 for supstep in range(schedule.supersteps_used + 1)}
    for node in schedule.T:
        nodesInSupstep[schedule.when[node]] += 1

    startIdx = schedule.supersteps_used - 1
    while startIdx >= 0:
        nodeCount = nodesInSupstep[startIdx]
        endIdx = startIdx
        while endIdx > 0 and (nodeCount + nodesInSupstep[endIdx - 1]) * schedule.instance.p * schedule.instance.p * (
                startIdx - endIdx + 1) <= batch_size:
            endIdx -= 1
            nodeCount += nodesInSupstep[endIdx]

        print(f'Next batch: supersteps {endIdx}-{startIdx}')
        partial = PartialParams(schedule.GetNodesInSupersteps(endIdx, startIdx), endIdx, startIdx)
        new_schedule, _ = run_ILP(name, "Partial", schedule, timeLimit, partial)
        if new_schedule is not None and new_schedule.objective_value < schedule.objective_value:
            schedule = new_schedule
        startIdx = endIdx - 1
        print(f'next val: {schedule.compute_objective()}')

    if schedule is not None:
        if not schedule.validate():
            print("ERROR: ILPIter returned invalid schedule.")
        else:
            print(f'Schedule: {schedule.compute_objective()}')
        solDict[name].add('ILPiter', schedule)
        schedule.write_file(f'{instance_name}', f'..{per}schedules{per}ILPiter_{name}')


# apply ILPfull and/or ILPpart
def ILP_apply(instance_names, solDict, timeLimit, batch_size=40):
    for instance_name in instance_names:
        name = os.path.split(instance_name)[1]
        inst = solDict[name].instance
        best_schedule = solDict[name].GetBestSchedule()
        print(f"Approximate number of variables in full ILP representation: {len(inst.G.nodes()) * inst.p * inst.p * best_schedule.supersteps_used}")
        if len(inst.G.nodes()) * inst.p * inst.p * best_schedule.supersteps_used < 20000:
            ILP_full_single(instance_name, solDict, 60)
        if "ILP" not in solDict[name].schedules or solDict[name].schedules["ILP"].objective_value > solDict[name].lower_bound:
            ILP_iter_single(instance_name, solDict, 3, batch_size)


def ILP_iter_all(instance_names, solDict, timeLimit, batch_size=40):
    for instance_name in instance_names:
        ILP_iter_single(instance_name, solDict, timeLimit, batch_size)


def ILP_CS_all(instance_names, solDict, timeLimit):
    for instance_name in instance_names:
        ILP_CS_single(instance_name, solDict, timeLimit)


# ILPinit
def ILP_init(instance_names, solDict, timeLimit, batch_size=2000):
    schedules_dir = f'..{per}schedules'
    for instance_name in instance_names:
        name = os.path.split(instance_name)[1]

        instance = solDict[name].schedules["Trivial"].instance
        list_init = {node: -1 for node in instance.G.nodes()}
        schedule = Schedule(instance, list_init, list_init, None)
        schedule.communication_schedule = {}
        top_order = instance.get_topological_order()
        superstep_per_batch = 3
        nr_of_nodes = ceil(batch_size / instance.p / instance.p / superstep_per_batch)
        stepIdx = 0
        for startIdx in range(0, len(top_order), nr_of_nodes):
            nodes = top_order[startIdx:min(startIdx+nr_of_nodes, len(top_order))]
            partial = PartialParams(nodes, stepIdx, stepIdx+superstep_per_batch-1)
            schedule_old = copy.deepcopy(schedule)
            schedule_partial = copy.deepcopy(schedule)
            schedule, _ = run_ILP(name, "Partial", schedule_partial, timeLimit, partial)
            if schedule is None:
                print("No schedule returned for subproblem; replacing it with trivial subschedule instead.")
                schedule = schedule_old
                for node in nodes:
                    schedule.where[node] = 0
                    schedule.when[node] = stepIdx
                    if startIdx+nr_of_nodes >= len(top_order):
                        schedule.invalid = False
            stepIdx += superstep_per_batch

        if schedule is not None:
            print(schedule.compute_objective())
            schedule.set_lazy_schedule()
            solDict[name].add('ILPInit', schedule)
            print(f'ILP init val: {schedule.compute_objective()}')
            schedule.write_file(f'{instance_name}', f'{schedules_dir}{per}ILPInit_{name}')

