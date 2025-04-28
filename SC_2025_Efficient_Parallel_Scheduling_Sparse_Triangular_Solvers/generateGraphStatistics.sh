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

gen_graph_stats() {
    if [ ! -f ${2} ]; then
        ${OSP_DIR}/build/examples/graph_analyser ${1} ${2} || exit 1
    fi
    echo "Finished ${2}."
    echo " "
}

OSP_DIR="$(pwd)/one-stop-parallel"

EVAL_DIR="$(pwd)/sc25-evaluation"

SUITESPARSE_DIR="$(pwd)/data_set/suite_sparse"
METIS_DIR="$(pwd)/data_set/metis"
ICHOL_DIR="$(pwd)/data_set/ichol"
ERDOS_RENYI_DIR="$(pwd)/data_set/erdos_renyi"
NARROW_BANDWIDTH_DIR="$(pwd)/data_set/narrow_bandwidth"

echo "This script generates all the graph statistics for the data analysis."

cd ${EVAL_DIR} || exit 1
gen_graph_stats ${NARROW_BANDWIDTH_DIR} graph_band_stats.txt
gen_graph_stats ${ERDOS_RENYI_DIR} graph_erdos_renyi_stats.txt
gen_graph_stats ${SUITESPARSE_DIR} graph_florida_stats.txt
gen_graph_stats ${ICHOL_DIR} graph_iChol_stats.txt
gen_graph_stats ${METIS_DIR} graph_hdagg_metis.txt

exit 0
