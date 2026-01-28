#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from ilp_parameters import instanceParameter, OBJ
from ilp import coptSolver 


def apply_params(p: instanceParameter, cfg: dict) -> None:
    # Core
    p.num_partitions = int(cfg["num_partitions"])
    p.memory_bound = float(cfg["memory_bound"])
    p.layer_threshold = int(cfg["layer_threshold"])
    p.per_layer_top_q = float(cfg.get("per_layer_top_q", 0.98))

    # model knobs
    p.cost_mla = float(cfg.get("cost_mla", 1.0))
    p.max_mlas_in_partition = int(cfg.get("max_mlas_in_partition", 10**9))
    p.lb_mla = int(cfg.get("lb_mla", 2))
    p.ub_mla = int(cfg.get("ub_mla", 4))
    p.load_threshold = float(cfg.get("load_threshold", 0.5))
    p.eps = float(cfg.get("eps", 1e-6))
    p.next_mla = bool(cfg.get("next_mla", False))

    # Time (solver multiplies by 3600)
    p.time_limit = float(cfg.get("time_limit_hours", 1.0))

    # Objective
    obj = str(cfg.get("objective", "LAMBDA")).upper()
    if obj == "LAMBDA":
        p.objective = OBJ.LAMBDA
    elif obj in {"FEAS", "FEASIBILITY", "NONE"}:
        # depends on enum; keep safe fallback
        p.objective = getattr(OBJ, "NONE", OBJ.LAMBDA)
    else:
        raise ValueError(f"Unknown objective: {obj}")

    # Optional warm start
    init_file = str(cfg.get("initial_solution_file", "")).strip()
    if init_file:
        p.initial_solution_file = init_file


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="pipeline.yaml")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        y = yaml.safe_load(f) or {}

    cfg = y.get("ilp")
    if not isinstance(cfg, dict):
        raise SystemExit("YAML missing ilp section")

    force = bool(cfg.get("force", False))
    pairs_pattern = str(cfg["pairs_pattern"])

    solutions_dir = Path(str(cfg.get("solutions_dir", "solutions")))
    solutions_dir.mkdir(parents=True, exist_ok=True)

    done_marker = solutions_dir / ".ilp_done"
    if done_marker.exists() and not force:
        print("[ilp] Done marker exists; skipping.")
        return 0

    p = instanceParameter()
    apply_params(p, cfg)

    base_file_name = str(cfg.get("base_file_name", str(solutions_dir / "partitioning_experts_pairs")))
    solver = coptSolver(p, base_file_name)

    transitions_folder = str(cfg.get("transitions_folder", "")).strip() or None

    solver.solve(pairs_pattern, transitions_folder)
    solver.write_solution_partition_to_file()

    done_marker.write_text("ok\n", encoding="utf-8")
    print(f"[ilp] Finished. Wrote done marker: {done_marker}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
