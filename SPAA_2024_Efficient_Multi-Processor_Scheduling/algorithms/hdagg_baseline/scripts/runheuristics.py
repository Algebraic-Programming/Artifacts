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
from schedule import read_schedule, CreateSolutionDict, Schedule
import instance
import util
from time import perf_counter
from sys import platform

per = "/" if platform == "linux" else "\\"


def run_heuristics(instance_names, solDict, settings):
    schedules_dirname = f'..{per}schedules'

    copy_instances(instance_names)
    prefix = "./" if platform == "linux" else ""

    os.chdir('..')
    if not os.path.isdir('temp_hdagg'):
        os.mkdir('temp_hdagg')
    os.chdir('temp_hdagg')
    copy = "cp" if platform == "linux" else "copy"
    os.system(f'{copy} ..{per}algos{per}lbc_demo {prefix}lbc_demo')

    for instance_path in instance_names:
        instance_name = os.path.split(instance_path)[1]
        mtx_name = instance_name + ".mtx"
        sched_name = instance_name + ".sched"
        solDict[instance_name].instance.write_DAG_mtx(mtx_name)
        try:
            os.system(f"{prefix}lbc_demo {prefix}{mtx_name} {solDict[instance_name].instance.p} > {sched_name}")
            lines = instance.GetCleanedFile(sched_name)
            processors = {entry[0]: entry[1] for entry in lines}
            supersteps = {entry[0]: entry[2] for entry in lines}
            solDict[instance_name].add('HDagg', Schedule(solDict[instance_name].instance, processors, supersteps, None))
        except Exception as e:
            print(f"Something went wrong with HDagg: {e}")

    os.chdir('..')
    os.system(f'rm -r temp_hdagg')
    os.chdir('scripts')



def copy_instances(instance_names):
    # copy instances
    copy = "cp" if platform == "linux" else "copy"
    for instance_path in instance_names:
        inst = os.path.split(instance_path)[1]
        os.system(f'{copy} {instance_path} ..{per}schedules{per}{inst}')


