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

from util import read_config, GetFileNames, CreateMergedFiles, RemoveMergedFiles
from runheuristics import run_heuristics, run_contraction, run_refinement
from runILP import ILP_apply, ILP_CS_all
from schedule import CreateSolutionDict
from os import path, remove


def Main_Run():
    config = read_config("config.txt")
    if config.file_mode == 'Combine':
        CreateMergedFiles(config)
    instances = GetFileNames(config)
    schedules = CreateSolutionDict(instances, config.settings)
    run_algorithms(instances, schedules, config.settings)
    if config.file_mode == 'Combine':
        RemoveMergedFiles()

    if path.isfile("results.CSV"):
        remove("results.CSV")
    names = [path.split(instance_path)[1] for instance_path in instances]
    for inst in names:
        print(f'\n{inst}:')
        schedules[inst].Print()
        schedules[inst].PrintCSV(inst, "results.CSV")


def run_algorithms(instances, schedules, settings):
    instance_names = instances
    solutions = schedules
    if settings.Contract:
        instance_names = run_contraction(instances, schedules, settings)
        solutions = CreateSolutionDict(instance_names, settings)
    run_heuristics(instance_names, solutions, settings)
    ILP_apply(instance_names, solutions, settings.TimeLimit, 4000)
    if settings.Contract:
        run_refinement(instances, solutions, schedules, settings)
    ILP_CS_all(instances, schedules, settings.TimeLimit)


if __name__ == '__main__':
    Main_Run()




