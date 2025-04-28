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
erdos_renyi
suite_sparse
ichol
metis
narrow_bandwidth
)

export OMP_NUM_THREADS=22
export OMP_PROC_BIND=close
export OMP_PLACES=cores

for arg in "${args[@]}"
do
    for f in data_set/$arg/*.mtx; do
        echo "Running $f"
        ./aggregation/build/example/Hdagg_SpTRSV $f 22
    done

    mv SpTrSv_Final_20.csv sc25-evaluation/SpTrSV_Data/SC_paper/SpTrSv_Final_20_$arg.csv

done

cp sc25-evaluation/SpTrSV_Data/SC_paper/SpTrSv_Final_20_suite_sparse.csv sc25-evaluation/SpTrSV_Data/plot_22_x86_SC/SpTrSv_Final_20_suite_sparse.csv