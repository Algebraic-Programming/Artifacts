#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Optional

import yaml
import coptpy as cp
from coptpy import COPT

from ilp_parameters import instanceParameter, OBJ
from hypergraph import Hypergraph
from hypergraph import MEM_MLA_IN_GB


class coptSolver(cp.CallbackBase):
    """
    Minimal solver.
    """

    def __init__(self, parameters: instanceParameter, base_file_name: str):
        super().__init__()
        self.parameters = parameters
        self.env = cp.Envr()

        self.model = None
        self.h: Optional[Hypergraph] = None

        # vars
        self.node_vars = None
        self.edge_vars = None
        self.t = None
        self.c = None
        self.z = None
        self.p = None

        self.num_edges = 0
        self.num_nodes = 0
        self.num_mla = 0

        self.base_file_name = base_file_name
        self.counter = 0

    # ---------- callback ----------
    def callback(self):
        if self.where() == COPT.CBCONTEXT_MIPSOL:
            self.write_solution_partition_to_file()
            self.counter += 1

    def write_solution_partition_to_file(self):
        if self.model is None:
            return

        with_mla = "with_next_mla" if self.parameters.next_mla else "no_next_mla"
        file_name = (
            f"{self.base_file_name}_"
            f"{self.parameters.num_partitions}_"
            f"{self.parameters.memory_bound}_"
            f"{self.parameters.layer_threshold}_"
            f"{with_mla}_"
            f"{self.counter}.txt"
        )
        Path(file_name).parent.mkdir(parents=True, exist_ok=True)

        with open(file_name, "w", encoding="utf-8") as fp:
            print(f"% Obj {self.model.objval}, # partitions {self.parameters.num_partitions}", file=fp)

            # Node -> partition assignment
            for i in range(self.num_nodes):
                for j in range(self.parameters.num_partitions):
                    if self.node_vars[i, j].x > 0.5:
                        print(f"{i} {j}", file=fp)

            layer_iter = range(self.parameters.layer_threshold)

            # MLA assignment per partition
            print("% --- MLA assignment per partition ---", file=fp)
            for j in range(self.parameters.num_partitions):
                layers_in_j = [l + 1 for l in layer_iter if self.t[l, j].x > 0.5]
                count = len(layers_in_j)
                layer_str = " ".join(map(str, layers_in_j)) if layers_in_j else "-"
                print(f"% partition {j}: #MLAs {count}; layers {layer_str}", file=fp)

            # Total MLAs per layer
            print("% --- Total MLAs per layer ---", file=fp)
            for l in layer_iter:
                tot = sum(1 for j in range(self.parameters.num_partitions) if self.t[l, j].x > 0.5)
                print(f"% layer {l + 1}: {tot}", file=fp)

    def set_initial_solution(self):
        if self.model is None:
            return

        sol_path = getattr(self.parameters, "initial_solution_file", "")
        if not sol_path or not os.path.exists(sol_path):
            return
        if self.h is None:
            return

        node_to_part = {}
        mla_layers_by_part = {}

        parsing_nodes = True
        mla_header = re.compile(r"%\s*---\s*MLA assignment per partition")
        mla_line = re.compile(r"%\s*partition\s+(\d+):.*layers\s+(.*)$")

        with open(sol_path, "r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue

                if line.startswith("%"):
                    if mla_header.search(line):
                        parsing_nodes = False
                    elif not parsing_nodes:
                        m = mla_line.match(line)
                        if m:
                            j = int(m.group(1))
                            rest = m.group(2).strip()
                            if rest in {"", "-"}:
                                mla_layers_by_part[j] = []
                            else:
                                mla_layers_by_part[j] = [int(tok) for tok in rest.split()]
                    continue

                if parsing_nodes:
                    try:
                        i_str, j_str = line.split()
                        node_to_part[int(i_str)] = int(j_str)
                    except Exception:
                        pass

        vars_list = []
        vals_list = []

        P = self.parameters.num_partitions
        L = self.num_mla

        # x[i,j]
        for i in range(self.num_nodes):
            assigned_j = node_to_part.get(i, None)
            for j in range(P):
                vars_list.append(self.node_vars[i, j])
                vals_list.append(1 if assigned_j == j else 0)

        # t[l,j] (file uses 1-based layers)
        t_one = set()
        for j, layers in mla_layers_by_part.items():
            for layer1 in layers:
                if 1 <= layer1 <= L:
                    t_one.add((layer1 - 1, j))

        for l in range(L):
            for j in range(P):
                vars_list.append(self.t[l, j])
                vals_list.append(1 if (l, j) in t_one else 0)

        self.model.setMipStart(vars_list, vals_list)
        self.model.loadMipStart()
        try:
            self.model.setParam(COPT.Param.MipStartMode, 2)
        except Exception:
            pass

    # ---------- ILP ----------
    def setup_moe_model_mla_v2(self, G: Hypergraph, para: instanceParameter):
        self.model = self.env.createModel("hgraph_part_moe_mla_v2")

        next_mla_constraints = para.next_mla

        # constants
        self.h = G
        self.num_nodes = G.num_nodes
        self.num_edges = G.num_edges
        self.num_mla = para.layer_threshold

        # variables
        self.node_vars = self.model.addVars(self.num_nodes, para.num_partitions, vtype=COPT.BINARY, nameprefix="x")
        self.edge_vars = self.model.addVars(self.num_edges, para.num_partitions, vtype=COPT.BINARY, nameprefix="y")
        self.t = self.model.addVars(self.num_mla, para.num_partitions, vtype=COPT.BINARY, nameprefix="t")
        self.z = self.model.addVars(self.num_mla, para.num_partitions, vtype=COPT.CONTINUOUS, lb=0.0, ub=1.0, nameprefix="z")
        
        self.c = self.model.addVars(self.num_mla, self.num_edges, para.num_partitions, vtype=COPT.BINARY, nameprefix="c")
        self.p = self.model.addVars(self.num_mla, self.num_edges, para.num_partitions, vtype=COPT.CONTINUOUS, lb=0.0, ub=1.0, nameprefix="p")

        # 1) partition constraint: sum_j x[i,j] = 1
        self.model.addConstrs(
            cp.quicksum(self.node_vars[i, j] for j in range(para.num_partitions)) == 1
            for i in range(self.num_nodes)
        )

        # 2) coupling of node and edge variables: y[e,j] >= x[h,j] for h in edge e
        self.model.addConstrs(
            self.edge_vars[e, j] >= self.node_vars[h, j]
            for j in range(para.num_partitions)
            for e in range(self.num_edges)
            for h in G.edges[e]
        )

        # 3) MLA-edge coupling + penalty scaffolding
        for l in range(self.num_mla):
            for e in range(self.num_edges):
                if l == G.edge_to_layer(e) - 1:
                    self.model.addConstrs(self.t[l, j] - self.edge_vars[e, j] <= self.c[l, e, j] for j in range(para.num_partitions))
                    self.model.addConstrs(self.p[l, e, j] + (1 - self.c[l, e, j]) >= self.z[l, j] for j in range(para.num_partitions))
                    if l + 1 < self.num_mla and next_mla_constraints:
                        self.model.addConstrs(self.t[l + 1, j] - self.edge_vars[e, j] <= self.c[l + 1, e, j] for j in range(para.num_partitions))
                        self.model.addConstrs(self.p[l + 1, e, j] + (1 - self.c[l + 1, e, j]) >= self.z[l + 1, j] for j in range(para.num_partitions))

        # 4) memory constraints
        self.model.addConstrs(
            cp.quicksum(G.node_mem_weights[i] * self.node_vars[i, j] for i in range(self.num_nodes))
            <= para.memory_bound - MEM_MLA_IN_GB * cp.quicksum(self.t[l, j] for l in range(self.num_mla))
            for j in range(para.num_partitions)
        )

        # 5) MLA replication bounds + load indicators
        for l in range(self.num_mla):
            expr = cp.quicksum(self.t[l, j] for j in range(para.num_partitions))
            self.model.addBoundConstr(expr, lb=para.lb_mla, ub=para.ub_mla)
            for j in range(para.num_partitions):
                self.model.addGenConstrIndicator(self.t[l, j], 1, self.z[l, j] >= para.load_threshold)
                self.model.addGenConstrIndicator(self.t[l, j], 0, self.z[l, j] <= para.eps)

        self.model.addConstrs(
            cp.quicksum(self.z[l, j] for j in range(para.num_partitions)) == 1.0
            for l in range(self.num_mla)
        )

        self.model.addConstrs(
            cp.quicksum(self.t[l, j] for l in range(self.num_mla)) <= para.max_mlas_in_partition
            for j in range(para.num_partitions)
        )

        # objective
        if para.objective == OBJ.LAMBDA:
            self.model.setObjective(
                cp.quicksum( G.edge_weights[e] * ( cp.quicksum(self.edge_vars[e, j] for j in range(para.num_partitions)) - 1) for e in range(self.num_edges))
                + 
                cp.quicksum( para.cost_mla
                            * ((self.p[l, e, j] if l >= 0 else 0.0) + (self.p[l + 1, e, j] if (l + 1 < self.num_mla and next_mla_constraints) else 0.0))
                            for j in range(para.num_partitions) for e in range(self.num_edges) for l in [G.edge_to_layer(e) - 1])                
                ,
                COPT.MINIMIZE,
            )
        else:
            # only LAMBDA supported here
            raise ValueError("This minimal ILP script supports only objective=OBJ.LAMBDA")

        return self.model

    def solve(self, pairs_pattern: str, folder_transitions: str | None = None):
        print(f"-- Solving (v2-only model) for pairs={pairs_pattern!r}")

        # Build graph (pairs-based)
        self.h = Hypergraph.from_multilayer_csvs_pairs(
            pairs_pattern,
            self.parameters.layer_threshold,
            self.parameters.per_layer_top_q,
            debug=True,
        )

        print(f"# edges {self.h.num_edges}, # nodes {self.h.num_nodes}")

        model = self.setup_moe_model_mla_v2(self.h, self.parameters)
        self.set_initial_solution()

        # solver params
        model.setParam(COPT.Param.TimeLimit, self.parameters.time_limit * 3600)
        model.setParam(COPT.Param.StrongBranching, 1)
        model.setParam(COPT.Param.SubMipHeurLevel, 0)
        model.setParam(COPT.Param.Presolve, 1)
        model.setParam(COPT.Param.CutLevel, 0)
        model.setParam(COPT.Param.TreeCutLevel, 2)
        model.setParam(COPT.Param.DivingHeurLevel, 2)
        model.setParam(COPT.Param.Threads, 192)
        model.setParam(COPT.Param.RelGap, 0.01)

        self.model.setCallback(self, COPT.CBCONTEXT_INCUMBENT)
        model.solve()


def apply_params(p: instanceParameter, cfg: dict) -> None:
    def get(k, default=None):
        return cfg[k] if k in cfg else default

    if "num_partitions" in cfg:
        p.num_partitions = int(cfg["num_partitions"])
    if "memory_bound" in cfg:
        p.memory_bound = float(cfg["memory_bound"])
    if "time_limit" in cfg:
        p.time_limit = float(cfg["time_limit"])  # hours
    if "layer_threshold" in cfg:
        p.layer_threshold = int(cfg["layer_threshold"])
    if "per_layer_top_q" in cfg:
        p.per_layer_top_q = float(cfg["per_layer_top_q"])

    if "cost_mla" in cfg:
        p.cost_mla = float(cfg["cost_mla"])
    if "eps" in cfg:
        p.eps = float(cfg["eps"])
    if "load_threshold" in cfg:
        p.load_threshold = float(cfg["load_threshold"])

    if "lb_mla" in cfg:
        p.lb_mla = int(cfg["lb_mla"])
    if "ub_mla" in cfg:
        p.ub_mla = int(cfg["ub_mla"])
    if "max_mlas_in_partition" in cfg:
        p.max_mlas_in_partition = int(cfg["max_mlas_in_partition"])
    if "next_mla" in cfg:
        p.next_mla = bool(cfg["next_mla"])

    if "initial_solution_file" in cfg:
        p.initial_solution_file = str(cfg["initial_solution_file"])

    # objective must be LAMBDA
    p.objective = OBJ.LAMBDA


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="pipeline.yaml")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        y = yaml.safe_load(f) or {}

    cfg = y.get("ilp")
    if not isinstance(cfg, dict):
        raise SystemExit("YAML missing ilp section")

    pairs_pattern = str(cfg["pairs_pattern"])
    base_file_name = str(cfg.get("base_file_name", "solutions/partitioning_experts_pairs"))
    solutions_dir = Path(cfg.get("solutions_dir", "solutions"))
    solutions_dir.mkdir(parents=True, exist_ok=True)

    force = bool(cfg.get("force", False))
    done_marker = solutions_dir / ".ilp_done"
    if done_marker.exists() and not force:
        print("[ilp] Done marker exists; skipping.")
        return 0

    p = instanceParameter()
    apply_params(p, cfg)

    solver = coptSolver(p, base_file_name)
    transitions_folder = str(cfg.get("transitions_folder", "")).strip() or None
    solver.solve(pairs_pattern, transitions_folder)
    solver.write_solution_partition_to_file()

    done_marker.write_text("ok\n", encoding="utf-8")
    print(f"[ilp] Finished. Wrote done marker: {done_marker}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
