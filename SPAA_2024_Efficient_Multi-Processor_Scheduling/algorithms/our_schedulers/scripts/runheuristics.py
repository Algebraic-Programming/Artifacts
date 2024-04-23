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

import os
from schedule import read_schedule, CreateSolutionDict, get_uncontracted_cost, get_uncontracted_schedule
import instance
import util
from time import perf_counter
from sys import platform
from runILP import ILP_init

per = "/" if platform == "linux" else "\\"


def run_heuristics(instance_names, solDict, settings):
    schedules_dirname = f'..{per}schedules'

    algo_names = get_algo_names(settings)

    if not settings.Contract:
        compile_algos(settings)
        copy_instances(instance_names)
        # in case of contraction, these were already done before

    # run baselines and initialization heuristics
    os.chdir(schedules_dirname)
    prefix = "./" if platform == "linux" else ""
    for instance_path in instance_names:
        instance = os.path.split(instance_path)[1]
        print(len(solDict[instance].instance.G.nodes()))
        for algo_name in algo_names:
            print(f"Running {algo_name} (plus HC) on {instance}")
            start = perf_counter()
            if algo_name == 'simple_schedulers':
                os.system(f"{prefix}{algo_name} -cilk -ETF -BL-EST -noHC -input {instance}")
                os.system(f"{prefix}{algo_name} -BSPg -timeLimit {settings.TimeLimit * 60} -input {instance}")

                solDict[instance].add('cilk', read_schedule(f'{schedules_dirname}{per}cilk_{instance}'))

                solDict[instance].add('BSPg', read_schedule(f'{schedules_dirname}{per}BSPg_{instance}'))
                solDict[instance].add('BSPgHCCS', read_schedule(f'{schedules_dirname}{per}BSPgHCCS_{instance}'))
                
                solDict[instance].add('ETF', read_schedule(f'{schedules_dirname}{per}ETF_{instance}'))
                solDict[instance].add('BL-EST', read_schedule(f'{schedules_dirname}{per}BL-EST_{instance}'))

            elif algo_name == 'second_heuristic_weights_cluster':
                try:
                    os.system(f"{prefix}{algo_name} {instance}")
                    solDict[instance].add('source3_weights_cluster',
                                          read_schedule(f'{schedules_dirname}{per}source3_weights_cluster_{instance}'))
                    # run HC on source3
                    os.system(
                        f"{prefix}simple_schedulers -onlyHC -timeLimit {settings.TimeLimit * 60} -input source3_weights_cluster_{instance}")
                    solDict[instance].add('source3_weights_clusterHCCS',
                                          read_schedule(f'{schedules_dirname}{per}source3_weights_clusterHCCS_{instance}'))
                except Exception as e:
                    print(f"Something went wrong with the clustering algorithm: {e}")

            print(f"    Took {round(perf_counter() - start, 2)} seconds")

    # run ILP initialization algorithm
    for instance_path in instance_names:
        name = os.path.split(instance_path)[1]
        if solDict[name].instance.p <= 4:
            try:
                ILP_init([instance_path], solDict, 2)
                os.system(f"{prefix}simple_schedulers -onlyHC -timeLimit {settings.TimeLimit * 60} -input ILPInit_{name}")
                solDict[name].add('ILPInitHCCS', read_schedule(f'{schedules_dirname}{per}ILPInitHCCS_{name}'))
            except Exception as e:
                print(f"Error during running ILPInit + HC: {e}")

    # delete executables and files
    if not settings.Contract:
        delete_files(instance_names, settings)
        # in case of contraction, this happens after refinement


# for multilevel algorithm
def run_contraction(instance_names, solDict, settings):
    compile_algos(settings)
    copy_instances(instance_names)
    os.chdir(f'..{per}schedules')
    prefix = "./" if platform == "linux" else ""
    contracted_instance_names = []

    # contract instances
    for instance_path in instance_names:
        instance_name = os.path.split(instance_path)[1]
        new_n = instance.compute_contraction_rate(solDict[instance_name].instance)
        os.system(f"{prefix}simple_schedulers -contract -shrinkTo {new_n} -input {instance_name}")
        contracted_instance_names.append(f'..{per}schedules{per}coarse_{instance_name}')
        print(f"Problem {instance_name} contracted to {new_n} nodes.")
    os.chdir(f'..{per}scripts')
    return contracted_instance_names


def run_refinement(instance_names, coarse_schedules, full_schedules, settings):
    os.chdir(f'..{per}schedules')
    prefix = "./" if platform == "linux" else ""
    # refine instances
    for instance_path in instance_names:
        name = os.path.split(instance_path)[1]
        contracted_name = f'contraction_{name}'
        #best_schedule = coarse_schedules["coarse_"+name].GetBestSchedule()
        best_schedule_name = coarse_schedules["coarse_"+name].GetBestName() +"_coarse_"+name
        print(best_schedule_name)
        HCsteps = 100
        os.system(f"{prefix}simple_schedulers -refine -input {name} -scheduleFile {best_schedule_name} -contractFile {contracted_name} -HCsteps {HCsteps}")
        print(f"{prefix}simple_schedulers -refine -input {name} -scheduleFile {best_schedule_name} -contractFile {contracted_name} -HCsteps {HCsteps}")
        full_schedules[name].add('refined', read_schedule(f'..{per}schedules{per}REF_{name}'))
        print(f"Problem {name} has been refined. Output schedule written to REF_{name}.")
        for algo, sched in coarse_schedules["coarse_"+name].schedules.items():
            if algo != "Trivial":
                full_schedules[name].add(algo, get_uncontracted_schedule(full_schedules[name].instance, coarse_schedules["coarse_"+name].schedules[algo], f'contraction_{name}'))
    os.chdir(f'..{per}scripts')
    delete_files(instance_names, settings)


# auxiliary functions
def compile_algos(settings):
    algo_names = get_algo_names(settings)
    # compile C++ sources
    os.chdir(f'..{per}algos')
    for algo_name in algo_names:
        os.system(f"g++ -o ..{per}schedules{per}{algo_name} {algo_name}.cpp")
    os.chdir(f'..{per}scripts')


def copy_instances(instance_names):
    # copy instances
    copy = "cp" if platform == "linux" else "copy"
    for instance_path in instance_names:
        inst = os.path.split(instance_path)[1]
        os.system(f'{copy} {instance_path} ..{per}schedules{per}{inst}')


def delete_files(instance_names, settings):
    os.chdir(f'..{per}schedules')
    prefix = "./" if platform == "linux" else ""
    suffix = "" if platform == "linux" else ".exe"
    algo_names = get_algo_names(settings)
    # delete C++ executables
    for algo_name in algo_names:
        os.remove(f"{prefix}{algo_name}{suffix}")
    # delete instance file copies
    for instance_path in instance_names:
        os.remove(f"{prefix}{os.path.split(instance_path)[1]}")
    os.chdir(f'..{per}scripts')


def get_algo_names(settings):
    algo_names = ['simple_schedulers']
    if settings.NUMA:
        algo_names.append('second_heuristic_weights_cluster')
    return algo_names
