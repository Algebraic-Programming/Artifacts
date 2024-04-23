LICENSE:
--------

Copyright 2024 Huawei Technologies Co., Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

SUMMARY:
--------

This directory contains efficient scheduling algorithms in the BSP model for a computation expressed in a hyperDAG format. It runs several initialization heuristics, then multiple hill climbing methods, and finally an ILP solver to obtain a schedule of small cost. It can consider different work/communication weights, different processor types, and NUMA effects between the processors.

Some of the algorithms only return a base BSP schedule that assigns each node to a given processor and superstep, while others also return a communication schedule (CS), which also describes when (i.e. in which superstep) we should send concrete values between given processors.

The program can be run by starting the file main.py in the ./scripts folder: just navigate to this folder, and type "python ./main.py" (in Linux). It is intended for Python versions 3.7-3.10.

As for the rest of the files and directories: the rest of the scripts controlling execution and the ILP formulation are in the ./scripts folder. The ./algo folder contains the C++ implementation of the initialization heuristics and the hill climbing algorithms. The output schedules are created in the ./schedules folder. The ./instances, ./merge1, ./merge2 folders contain different kinds of sample problem inputs.

Disclaimer: the tool is a research prototype, and as such, it might occasionally contain unexpected bugs or runtime errors. It also assumes that it is run with a correct parametrization: the correctness of the program parameters, and the correct format of the input file are only checked to a limited extent.

CONFIG FILE:
------------

The config settings are obtained from the config.txt file in the ./scripts folder. If this file has an incorrect format or there is a read error, the default settings are used instead.

The config file can contain any number of comment lines, starting with the character "%".
The rest of the lines must contain a valid parameter assignment, but these can be in any order. The parameter assignment lines can either contain 2 words (<parameter> <value>), or 3 words with an equation sign in the middle (<parameter> = <value>).

The script can be run in 3 file modes: File, Dir or Combine. In case of File, the script is run on the given problem instance. In case of Dir, the script reads all the files in a given directory, attempting to interpret each as a problem instance. In case of Combine, the script takes two directories as an input, with the DAG descriptions in the first directory, and the machine parameters in the second one; the script combines and solves all the problem combinations from these input DAGs and machines. The default parametrization is "Dir = ../instances".
Some examples lines for parametrizations are: "Dir = ../test", "File = ../instances/example.txt" or "Combine = ../merge1 ../merge2"

The TimeLimit parameter is set to 3 by default. It needs to be a positive integer, and it specifies the time limit (in minutes) for the ILPcs and hill climbing algorithms that are run by the script.
Exmaple line: "TimeLimit = True"

The Contract parameter is set to False by default. If set to true, then the multilevel scheduling algorithm is applied (which is a somewhat preliminary implementation of this approach).
Exmaple line: "Contract = True"


OUTPUT SCHEDULES:
-----------------

The output schedules are all provided in the ./schedules folder. For any problem instance, some subset of the following output schedules is created:
- cilk (as a baseline)
- ETF and BL-EST (as baselines)
- BSPg, source3_weights_cluster and potentially ILPInit (initialization heuristics), also improved afterwards with hill climbing (+HC), and both standard and CS-specific hill climbing (+HCCS).
- ILP and/or ILPpart (here called ILPiter) and/or ILPcs (corresponding to a general ILP, the partial ILP formulation and a CS-specific ILP).

Besides this, a summary is also displayed both on the standard output and in a results.CSV file in the ./scripts folder, containing the cost of the found schedules (plus the number of supersteps and the best lower bound) for each instance and each algorithm.


FILE FORMATS:
-------------

The input problem instances are text files, and consist of a (i) DAG description and a (ii) machine description, in this order. In case of "Combine" mode, these two parts are separated into two distinct files. The output schedule files, on the other hand, first consists of the same problem instance description, which is then followed by the output schedule.
Besides the predefined file structure, the files may also contain any number of comment lines at any point; each of these must start with a "%" symbol. All lines are ended by the line feed character ('\n').

1. DAG description:

- The first line of the DAG description contains three integers 'H N M', separated by spaces. 'H' describes the number of hyperedges, 'N' describes the number of nodes, and 'M' describes the number of pins (i.e. sum of the total size of all hyperedges).

- This is then followed by 'M' distinct lines, each describing a specific pin of the hyperDAG. Each of these lines contain two integers 'E V', separated by a space. Integer 'E' lies between '0' and 'H-1' (inclusive) while 'V' lies between '0' and 'N-1'. The line indicates a pin connecting hyperedge 'E' to node 'V'. The first pin listed for each hyperedge describes the *source node* of the given hyperedge in the hyperDAG: in the original DAG representation, there is a directed edge from this source node to all the other nodes contained in this hyperedge.

- This is then followed by 'N' distinct lines, each describing the properties of a specific node. Each of these lines begins with the index of a node, i.e., an integer between '0' and 'N-1'; each of the 'N' lines must begin with a unique integer in this range. Each line can contain further integers, separated by spaces, describing properties of this node. The second number in the line describes the work weight of the node (the time required to compute it), while the third number describes the communication weight of the node (i.e. the size of the output, in case it needs to be sent to other processors). Both of these weights must be non-negative integers.

2. Machine description / problem parameters:

- The first line of the machine description contains three integers 'P G L'. 'P' describes the number of processors, 'G' describes the communication cost of a single unit of data, and 'L' describes the latency cost incurred by each new superstep.

- In case we have NUMA data in the file format (which is always the case in the current examples), the next P*P lines describe the NUMA cost between each specific pair of processors (note that these are only relative costs; for the concrete cost, each of them will be multiplied by 'G'). Each of these P*P lines must contain three integers 'P1 P2 C', indicating that the relative cost of sending data from processor 'P1' to processor 'P2' is 'C'. The first two integers 'P1' and 'P2' are processor indices between '0' and 'P-1'; for each valid ordered pair of indices between '0' and 'P-1', we must have a unique such line. 'C' is a non-negative integer that describes the cost of sending data from 'P1' to 'P2'.

3. Schedule description:

The output files begin with the problem description above (DAG and machine description) and then continue with the following schedule description. 

- It first begins by 'N' distinct lines, each describing a specific node. Each of these lines contain three integers 'V P1 S1', separated by spaces, indicating that in the solution, node 'V' is assigned to processor 'P1' and to superstep 'S1'. The indices 'V' must be unique integers between '0' and 'N-1' in each of these 'N' lines. 'P1' must be a processor index between '0' and 'P-1'. 'S1' is a non-negative integer, the index of a superstep.

- Optionally, the schedule description may continue with a communication schedule; this is e.g. contained in the output when ILP, ILPCS or HCCS is used, but not in BSPg, source3_weights_cluster and HC. The communication schedule begins with a line containing a single integer 'Q', describing the number of following lines. This is followed by 'Q' distinct lines, describing the communication schedule. Each of these 'Q' lines consist of four integers 'V P1 P2 S', separated by spaces, indicating that the output of node 'V' is transferred from processor 'P1' to processor 'P2' in the communication phase of superstep 'S'. 'V' is a node index between '0' and 'N-1', 'P1' and 'P2' are distinct processor indices between '0' and 'P-1', and 'S' is a non-negative integer indicating a superstep index.


