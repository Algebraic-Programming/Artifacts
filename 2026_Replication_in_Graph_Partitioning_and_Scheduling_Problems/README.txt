This repository is copyright by the Computing Systems Laboratory, Zurich Research Center, Huawei Technologies Switzerland AG.

Data and tools are distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the applicable licenses for the specific language governing permissions and limitations.

# Summary

The repository contains the algorithms and data corresponding to the paper "Replication in Graph Partitioning and Scheduling Problems", by Pál András Papp, Toni Böhnlein and Albert-Jan Yzelman, in 2026. The full version of the paper is available on arXiv. The role of this directory is to ensure reproducibility for our empirical results in the paper above.

# License

All tools in this directory are licensed under the Apache License, Version 2.0 (the "License"); you may not use the tools except in compliance with the License. You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0


# Algorithms

The implementation of the algorithms in the paper is available as part of the OneStopParallel scheduling toolbox (OSP), available at https://github.com/Algebraic-Programming/OneStopParallel. The experiments were conducted with the version of this codebase from 18th February 2026.

The partitioning ILPs are part of a separate partitioning module in OneStopParallel, whereas the BSP scheduling ILPs and heuristics are incorporated into the core scheduling module. The code contains a partitioning test suite for the different ILP algorithms (on a separate partitioning_test_suite branch at the time of writing), and the BSP scheduling heuristics are part of the BSP test suite. In the latter, the different variants of the advanced replication heuristic can be accessed via setting a method mask to the greedy replicating algorithm, where 1=batch replication, 2=superstep merging and 4=superstep replication.

For more details on the algorithms, we refer to the paper or the documentation of OneStopParallel.

# Datasets

For partitioning, there are four datasets.

- Two of them originate from the SuiteSparse sparse matrix collection available at https://sparse.tamu.edu/. OneStopParallel can directly read these inputs; the corresponding SpMV hypergraph model is determined by the file extension.

- The other two datasets originate from Mixture of Experts LLM serving, and were created as discussed in the paper, from the data available for the paper "Orders in Chaos: Enhancing Large-Scale MoE LLM Serving with Data Movement Forecasting", by Zhongkai Yu, Yue Guan, Zihao Yu, Chenyang Zhou, Zhengding Hu, Shuyi Pei, Yangwook Kang, Yufei Ding, and Po-An Tsai, available as arXiv preprint arXiv:2510.05497 (2025). We directly upload these hypergraphs here in the MoE_dataset subdirectory, since this benchmark is one of the new contributions of the paper. The hypergraphs here are in a .hmetis hypergraph format, and can also be directly read by OSP.

For scheduling, there are three datasets:

- The hdb dataset comes form our earlier paper "Efficient Multi-Processor Scheduling in Increasingly Realistic Models", published in the 36th ACM Symposium on Parallelism in Algorithms and Architectures, 2024. The corresponding computational DAGs are available at github.com/Algebraic-Programming/Artifacts/tree/master/SPAA_2024_Efficient_Multi-Processor_Scheduling , and can be read again with OSP.

- The PSDD dataset originates from the StarAI Circuit Model Zoo at https://github.com/UCLA-StarAI/Circuit-Model-Zoo/tree/master/psdds. These files are in a custom format, and hence need to be first converted to a format readable with OSP.

- The sptrsv dataset is again derived from sparse matrices from the SuiteSparse sparse matrix collection at https://sparse.tamu.edu/.

Some examples for machine parameters files for scheduling can also be found in the above mentioned artifact folder for the SPAA '24 paper.

# Result data

Finally, the result_data folder contains the raw data combined form the output files of the OSP test suites; our experimental results are based on this data.

For hypergraph partitioning, the results are contained in hypergraph_partitioning_results.csv. Each line here corresponds (one of the three) replicating ILP formulations run on a concrete instance with a concrete number of partitions and balance parameter. Note that the files also contain some duplicate entries; these typically stem from cases when the server or the ILP solver encountered a problem, so the solving process needed to be repeated. The ablation study data for hypergraph partitioning is contained in a similar format in a separate file, hypergraph_partitioning_ablation.csv.

For scheduling, the heuristics' results are contained in bsp_scheduling_heuristics_results.csv. Each line here corresponds to a concrete instance, parameter combination and algorithm; the latter is either the GreedyBspHC baseline, or our replicating heuristic with a concrete mask for the enabled methods.

For scheduling, the ILP results are contained in bsp_scheduling_ilp_results.csv. In contrast to previous files, this data contains the cost of both the replicating and the non-replicating solution (for a concrete computational DAG and parameter configuration) in the same row, to allow for a simpler comparison.