This repository is copyright by the Computing Systems Laboratory, Zurich Research Center, Huawei Technologies Switzerland AG.

Data and tools are distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the applicable licenses for the specific language governing permissions and limitations.

# Summary

The repository contains our DAG scheduling framework, which combines several algorithms to efficiently schedule arbitrary computational DAGs in the BSP model.

The algorithms and experiments are discussed in detail in our paper "Efficient Multi-Processor Scheduling in Increasingly Realistic Models", published in the 36th ACM Symposium on Parallelism in Algorithms and Architectures, 2024, or in the full version of the paper available on arXiv.

We note that the role of this directory is to ensure reproducibility for our empirical results in the paper above. For a significantly upgraded version of our scheduling algorithms, consider our scheduling tool available in a folder within https://github.com/Algebraic-Programming/OneStopParallel.

# License

All tools in this directory are licensed under the Apache License, Version 2.0 (the "License"); you may not use the tools except in compliance with the License. You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0


# Structure

## ./algorithms

This folder contains the source codes of our algorithms, as well as the necessary tools to run some baselines.

## ./database

This folder contains the computational DAG and machine parameter files used throughout our experiments.

## ./experiment_results

This folder contains the output data from our experiments, as welll as a simple script to process and aggregate these data files.

