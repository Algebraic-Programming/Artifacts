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

from dataclasses import dataclass
from itertools import groupby
from random import randint
from typing import Dict, Tuple
from collections import defaultdict
from math import ceil

import networkx as nx


@dataclass
class Instance:
    G: nx.classes.digraph.DiGraph
    L: float  # latency cost of superstep
    g: float  # scaling factor communication cost
    p: float  # number of processors
    numa_costs: Dict[Tuple[int, int], int] = None  # cost factor to communicate between p1 and p2
    recomputation: bool = False

    def __post_init__(self):
        # set attributes if they were not given
        if self.numa_costs is None:
            self.numa_costs = NUMA_default(self.p)
        if nx.get_node_attributes(self.G, 'compute_cost') == {}:
            compute_costs = {v: randint(1, 10) for v in self.G.nodes()}
            nx.set_node_attributes(self.G, compute_costs, name="compute_cost")
        if nx.get_node_attributes(self.G, 'communication_cost') == {}:
            communication_costs = {v: randint(1, 3) for v in self.G.nodes()}
            nx.set_node_attributes(self.G, communication_costs, name="communication_cost")

    def get_topological_order(self):
        out_edges = {node: [] for node in self.G.nodes()}
        in_degree = {node: 0 for node in self.G.nodes()}
        for source, target in self.G.edges():
            out_edges[source].append(target)
            in_degree[target] += 1

        top_order = []
        predecessor_count = {node: 0 for node in self.G.nodes()}
        queue = []

        for node in self.G.nodes():
            if in_degree[node] == 0:
                queue.append(node)

        while len(queue):
            node = queue[0]
            queue.pop(0)
            top_order.append(node)
            for succ in out_edges[node]:
                predecessor_count[succ] += 1
                if predecessor_count[succ] == in_degree[succ]:
                    queue.append(succ)

        return top_order
        

    def write_DAG_mtx(self, filename):
        out_edges = {node: [] for node in self.G.nodes()}
        for source, target in self.G.edges():
            out_edges[source].append(target)
        total = 0
        for node in self.G.nodes():
            out_edges[node].append(node)
            out_edges[node].sort()
            total += len(out_edges[node])
        with open(filename, 'w') as file:
            file.write('%%MatrixMarket matrix coordinate integer symmetric\n')
            file.write(f'{len(self.G.nodes())} {len(self.G.nodes())} {total}\n')
            for node in self.G.nodes():
                for succ in out_edges[node]:
                    file.write(f'{succ+1} {node+1} 1\n')
        print(f"DAG mtx written to {filename}")


def GetCleanedFile(filename):
    cleaned = []
    with open(filename, 'r') as file:
        for line in file:
            clean_line = line.strip().split()
            if line.startswith('%') or len(clean_line) == 0:
                continue
            else:
                cleaned.append([int(string) for string in clean_line])
    return cleaned


def GetInstanceFromCleanedFile(cleaned, NUMA=True):
    num_hyperedges, num_nodes, num_pins = cleaned[0]
    pins = cleaned[1:1 + num_pins]
    nodes = cleaned[num_pins + 1:num_pins + 1 + num_nodes]

    p, g, L = cleaned[num_pins + num_nodes + 1]

    numa_lines = p * p if NUMA else 0
    if NUMA:
        numa_data = cleaned[num_pins + num_nodes + 2:num_pins + num_nodes + 2 + numa_lines]
        numa_costs = {(p1, p2): entry for p1, p2, entry in numa_data}
    else:
        numa_costs = NUMA_default(p)

    node_info = [(index, {'compute_cost': compute_cost, 'communication_cost': communication_cost})
                 for index, compute_cost, communication_cost, in nodes]

    nets = [[entry[1] for entry in list(group)] for key, group in groupby(pins, lambda x: x[0])]

    G = nx.DiGraph()
    G.add_nodes_from(node_info)
    G.add_edges_from([(net[0], net[i]) for net in nets for i in range(1, len(net))])

    return Instance(G, L, g, p, numa_costs)


def read_instance(filename, NUMA=True):
    cleaned = GetCleanedFile(filename)
    return GetInstanceFromCleanedFile(cleaned, NUMA)


def NUMA_default(p):
    numa_costs = {(p1, p2): 1 for p1 in range(p) for p2 in range(p)}
    for p1 in range(p):
        numa_costs[(p1, p1)] = 0
    return numa_costs


def read_DAG(filename, num_entries):
    def clean(line):
        entries = line.strip().split()
        length = len(entries) == 2
        types = all([entry.isdigit() for entry in entries])
        return all([length, types])

    def read(filename, num_entries):
        counter = 0
        data = defaultdict(list)
        with open(filename, 'r') as outfile:
            for line in outfile.readlines():
                if clean(line) and counter <= num_entries:
                    counter += 1
                    h_id, node = line.strip().split()
                    data[h_id].append(node)
            data = {h_id: net for h_id, net in data.items() if len(net) > 1}  # remove singleton edges that can occur if reading is cut off
            return data

    data = read(filename, num_entries)
    G = nx.DiGraph()
    G.add_edges_from([(net[0], net[i]) for net in data.values() for i in range(1, len(net))])
    return G


