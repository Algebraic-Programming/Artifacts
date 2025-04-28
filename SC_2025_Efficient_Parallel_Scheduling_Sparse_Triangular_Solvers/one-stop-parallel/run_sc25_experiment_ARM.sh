#!/bin/bash

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

# @author Toni Boehnlein, Pal Andras Papp, Raphael S. Steiner

export OMP_NUM_THREADS=22
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Calling ./build/test_suite/test_suite_execution_SM --config sc25-configs/main_experiment_ARM.json"
    ./build/test_suite/test_suite_execution_SM --config sc25-configs/main_experiment_ARM.json

echo "Calling ./build/test_suite/test_suite_execution --config sc25-configs/serial_ARM.json"
    ./build/test_suite/test_suite_execution --config sc25-configs/serial_ARM.json




