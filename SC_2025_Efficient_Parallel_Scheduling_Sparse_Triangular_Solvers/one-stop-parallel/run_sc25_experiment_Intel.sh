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

args=(
sc25-configs/main_experiment_suite_sparse_x86.json
sc25-configs/main_experiment_erdos_renyi_x86.json
sc25-configs/main_experiment_ichol_x86.json
sc25-configs/main_experiment_metis_x86.json
sc25-configs/main_experiment_narrow_bandwidth_x86.json
sc25-configs/threads_experiments_x86.json
)

export OMP_NUM_THREADS=22
export OMP_PROC_BIND=close
export OMP_PLACES=cores

for arg in "${args[@]}"
do
    echo "Calling ./build/test_suite/test_suite_execution_SM --config $arg"
    ./build/test_suite/test_suite_execution_SM --config $arg
done


args_serial=(
sc25-configs/serial_suite_sparse_x86.json
sc25-configs/serial_erdos_renyi_x86.json
sc25-configs/serial_ichol_x86.json
sc25-configs/serial_metis_x86.json
sc25-configs/serial_narrow_bandwidth_x86.json
)

for arg in "${args_serial[@]}"
do
    echo "Calling ./build/test_suite/test_suite_execution --config $arg"
    ./build/test_suite/test_suite_execution --config $arg
done

echo "Calling ./build/test_suite/test_suite_execution_SM_no_perm --config sc25-configs/reorder_experiment_x86.json"
./build/test_suite/test_suite_execution_SM_no_perm --config sc25-configs/reorder_experiment_x86.json

cp ../sc25-evaluation/SpTrSV_Data/SC_paper/run_serial_suite_sparse.csv ../sc25-evaluation/SpTrSV_Data/plot_22_reorder/run_serial_suite_sparse.csv
cp ../sc25-evaluation/SpTrSV_Data/SC_paper/all_run_suite_sparse.csv ../sc25-evaluation/SpTrSV_Data/plot_22_reorder/all_run_suite_sparse.csv

cp ../sc25-evaluation/SpTrSV_Data/SC_paper/run_serial_suite_sparse.csv ../sc25-evaluation/SpTrSV_Data/plot_22_x86_SC/run_serial_suite_sparse.csv
cp ../sc25-evaluation/SpTrSV_Data/SC_paper/all_run_suite_sparse.csv ../sc25-evaluation/SpTrSV_Data/plot_22_x86_SC/all_run_suite_sparse.csv