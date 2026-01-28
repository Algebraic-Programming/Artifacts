from __future__ import annotations
from collections import Counter
from pathlib import Path
import csv
import os
from typing import List, Dict, FrozenSet, Tuple, Union
from utility_hyperedges import *

KeyType = Union[FrozenSet[int], Tuple[int, ...]]

from dotenv import load_dotenv

load_dotenv()

def count_by_layer_from_file(
    path: str | Path,
    k: int = 4,
    ignore_negatives: bool = True,
    as_tuples: bool = False,
    debug: bool = False,
) -> List[Dict[KeyType, int]]:
    """
    Single-pass reader, updating per-layer hyperedge counts on the fly.
    Assumes each layer line has exactly k experts.
    """
    if debug:
        set_debug(True)
        logger.debug(f"Starting count_by_layer_from_file(path={path}, k={k}, "
                     f"ignore_negatives={ignore_negatives}, as_tuples={as_tuples})")

    layers: List[Counter[KeyType]] = [Counter() for _ in range(N_LAYERS)]
    key_builder = (lambda seq: tuple(seq)) if as_tuples else (lambda seq: frozenset(seq))

    with open(path, "r", encoding="utf-8") as f:
        for line_no, raw in enumerate(f, start=1):
            sraw = raw.strip()

            if not sraw:
                if debug: logger.debug(f"[line {line_no}] blank -> skip")
                continue

            if DIVIDER_RE.match(sraw):
                if debug: logger.debug(f"[line {line_no}] divider -> skip")
                continue

            m = LAYER_LINE_RE.match(sraw)
            if not m:
                if debug: logger.debug(f"[line {line_no}] no match -> skip: {sraw!r}")
                continue

            layer_idx = int(m.group(1))  # 1-based index
            if not (1 <= layer_idx < N_LAYERS):
                if debug: logger.debug(f"[line {line_no}] layer {layer_idx} out of range -> skip")
                continue

            # parse experts
            try:
                # parse experts (handles "id" or "id:score")
                experts = []
                for x in m.group(2).split(","):
                    tok = x.strip()
                    if not tok:
                        continue
                    tok = tok.split(":", 1)[0].strip()  # keep left side if "id:score"
                    experts.append(int(tok))

                if ignore_negatives and -1 in experts:
                    if debug:
                        logger.debug(f"[line {line_no}] contains -1 -> skip (experts={experts})")
                    continue

                if len(experts) != k:
                    if debug:
                        logger.debug(f"[line {line_no}] expected k={k}, got {len(experts)} -> skip")
                    continue

                key = key_builder(sorted(experts))
                layers[layer_idx][key] += 1

                if debug:
                    logger.debug(
                        f"[line {line_no}] layer={layer_idx} key={tuple(sorted(experts))} "
                        f"-> freq={layers[layer_idx][key]}"
                    )
            except ValueError as e:
                if debug:
                    logger.debug(f"[line {line_no}] parse error: {e} in {m.group(2)!r} -> skip line")
                continue

            

    result = [dict(c) for c in layers]

    if debug:
        totals = sum(sum(d.values()) for d in result)
        logger.debug(f"Finished. Nonempty layers: {sum(1 for d in result if d)}; total increments: {totals}")
        set_debug(False)  # restore default

    return result

# ---------------- saving ----------------
def _serialize_key(key: KeyType) -> Tuple[int, ...]:
    """Make keys CSV friendly (sorted tuple of ints)."""
    if isinstance(key, tuple):
        return key
    return tuple(sorted(key))  # frozenset -> sorted tuple

def save_per_layer(
    per_layer: List[Dict[KeyType, int]],
    out_dir: str | Path,
    min_freq: int = 1,
    file_prefix: str = "layer",  # e.g., "layer_01.csv"
    zero_pad: int = 2,
    skip_empty: bool = True,
    debug: bool = False,
) -> None:
    """
    Save one CSV per layer in out_dir, named like layer_01.csv, layer_02.csv, ...
    """
    if debug:
        set_debug(True)
        logger.debug(f"Saving per-layer CSVs to {out_dir} (min_freq={min_freq})")

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    written = 0
    for idx, d in enumerate(per_layer, start=1):
        rows = [(_serialize_key(k), v) for k, v in d.items() if v >= min_freq]
        if not rows and skip_empty:
            if debug: logger.debug(f"[layer {idx}] empty after threshold -> skip")
            continue

        filename = f"{file_prefix}_{str(idx).zfill(zero_pad)}.csv"
        fp = out_path / filename

        if debug:
            logger.debug(f"[layer {idx}] writing {len(rows)} rows -> {fp}")

        with fp.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["experts", "frequency"])
            # sort by freq desc, then lexicographically
            for experts, freq in sorted(rows, key=lambda x: (-x[1], x[0])):
                w.writerow([";".join(map(str, experts)), freq])

        written += 1

    if debug:
        logger.debug(f"Done. Files written: {written}")
        set_debug(False)



if __name__ == "__main__":
    import argparse
    import yaml

    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="pipeline.yaml")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    c = cfg.get("extract_hyperedge", {})

    out_dir = Path(c["out_dir"])
    file_prefix = str(c.get("file_prefix", "layer"))
    zero_pad = int(c.get("zero_pad", 2))
    skip_empty = bool(c.get("skip_empty", True))
    force = bool(c.get("force", False))  # optional YAML switch

    # If skip_empty=True, we can't know exactly which layers will be written without running.
    # So we use a conservative "already done" check:
    # - If there is at least one matching layer CSV in out_dir, assume stage is done.
    # For a stricter check, set skip_empty: false in YAML.
    if not force:
        if skip_empty:
            existing = sorted(out_dir.glob(f"{file_prefix}_*.csv"))
            if existing:
                print(f"[extract_hyperedge] Found existing outputs in {out_dir} (e.g. {existing[0].name}); skipping.")
                raise SystemExit(0)
        else:
            # strict check: all layer_XX.csv for 1..N_LAYERS-1 exist
            expected = [
                out_dir / f"{file_prefix}_{str(i).zfill(zero_pad)}.csv"
                for i in range(1, N_LAYERS) 
            ]
            if all(p.exists() for p in expected):
                print(f"[extract_hyperedge] All expected outputs already exist in {out_dir}; skipping.")
                raise SystemExit(0)

    # ---- run stage normally ----
    path = c["input_path"]
    input_path = os.path.expandvars(path)
    if "$" in input_path:
        raise RuntimeError(f"Unresolved env var in input_path: {path}")

    if not os.path.exists(input_path):
        raise FileNotFoundError(input_path)
    per_layer = count_by_layer_from_file(
        input_path,
        k=int(c.get("k", K)),
        ignore_negatives=bool(c.get("ignore_negatives", True)),
        as_tuples=bool(c.get("as_tuples", False)),
        debug=bool(c.get("debug", False)),
    )

    save_per_layer(
        per_layer,
        out_dir=out_dir,
        min_freq=int(c.get("min_freq", MIN_FREQ)),
        file_prefix=file_prefix,
        zero_pad=zero_pad,
        skip_empty=skip_empty,
        debug=bool(c.get("debug", False)),
    )
