#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import itertools
import os
import re
from collections import Counter
from pathlib import Path

import yaml


def read_rows(file_path: str, experts_col: str, freq_col: str, encoding: str):
    with open(file_path, "r", newline="", encoding=encoding) as f:
        reader = csv.DictReader(f)
        for row in reader:
            experts_field = (row.get(experts_col) or "").strip()
            freq_field = (row.get(freq_col) or "").strip()
            if not experts_field or not freq_field:
                continue

            experts = []
            for s in experts_field.split(";"):
                try:
                    experts.append(int(s.strip()))
                except ValueError:
                    pass

            if len(experts) < 2:
                continue

            experts = sorted(set(experts))
            try:
                freq = int(freq_field)
            except ValueError:
                continue

            yield experts, freq


def process_file(
    file_path: str,
    output_dir: str,
    experts_col: str,
    freq_col: str,
    input_encoding: str,
    output_encoding: str,
    *,
    force: bool,
):
    os.makedirs(output_dir, exist_ok=True)

    base = os.path.basename(file_path)
    name, _ = os.path.splitext(base)
    out_file = os.path.join(output_dir, f"{name}_pairs.csv")

    if (not force) and os.path.exists(out_file):
        print(f"[extract_pairs] {out_file} exists; skipping.")
        return

    pair_counts = Counter()
    for experts, freq in read_rows(file_path, experts_col, freq_col, input_encoding):
        for a, b in itertools.combinations(experts, 2):
            pair_counts[(a, b)] += freq

    with open(out_file, "w", newline="", encoding=output_encoding) as out:
        writer = csv.writer(out)
        writer.writerow(["expert_a", "expert_b", "frequency"])
        for (a, b), total in sorted(pair_counts.items()):
            writer.writerow([a, b, total])

    print(f"Wrote {len(pair_counts)} pairs to {out_file}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="pipeline.yaml")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    c = cfg.get("extract_pairs")
    if not isinstance(c, dict):
        raise SystemExit("YAML missing extract_pairs section")

    input_pattern = str(c["input_pattern"])
    output_dir = str(c["output_dir"])
    valid_re = re.compile(str(c["valid_file_regex"]))

    experts_col = str(c.get("experts_col", "experts"))
    freq_col = str(c.get("frequency_col", "frequency"))
    input_encoding = str(c.get("input_encoding", "utf-8-sig"))
    output_encoding = str(c.get("output_encoding", "utf-8"))
    force = bool(c.get("force", False))

    files = sorted(glob.glob(input_pattern))
    files = [f for f in files if valid_re.match(os.path.basename(f))]
    if not files:
        print(f"No files matched pattern: {input_pattern}")
        return 1

    # skip whole step if ALL outputs already exist
    if not force:
        expected = [Path(output_dir) / (Path(f).stem + "_pairs.csv") for f in files]
        if expected and all(p.exists() for p in expected):
            print("[extract_pairs] All pair files already exist; skipping step.")
            return 0

    for fp in files:
        process_file(
            fp,
            output_dir,
            experts_col,
            freq_col,
            input_encoding,
            output_encoding,
            force=force,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
