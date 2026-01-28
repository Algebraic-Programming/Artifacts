#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, Tuple, List, Iterator, Set, Optional

import yaml
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dotenv import load_dotenv

load_dotenv()


# ----------------------- Streaming parser -----------------------

def stream_layers(path: str) -> Iterator[Tuple[int, int, List[int]]]:
    """
    Yield (chunk_id, layer_id, experts_list).
    Chunks are separated by lines consisting only of '%' characters.
    Ignores '-1' entries.
    """
    chunk_id = 0
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s[0] == '%' and all(c == '%' for c in s):
                chunk_id += 1
                continue

            i = s.find(':')
            if i < 0:
                continue
            lid_str = s[:i].strip()
            try:
                layer_id = int(lid_str)
            except ValueError:
                continue

            rest = s[i + 1:]
            exps: List[int] = []
            for tok in rest.split(','):
                tok = tok.strip()
                if not tok:
                    continue
                try:
                    v = int(tok)
                except ValueError:
                    continue
                if v != -1:
                    exps.append(v)

            yield (chunk_id, layer_id, exps)

# ----------------------- Dict accumulator -----------------------

from collections import defaultdict

def accumulate_layerpair_edges_dict(path: str):
    edges_by_pair: Dict[Tuple[int,int], Dict[Tuple[int,int], int]] = defaultdict(lambda: defaultdict(int))
    sources_by_pair: Dict[Tuple[int,int], Set[int]] = defaultdict(set)
    targets_by_pair: Dict[Tuple[int,int], Set[int]] = defaultdict(set)

    prev_chunk: Optional[int] = None
    prev_layer: Optional[int] = None
    prev_exps: Optional[List[int]] = None

    for chunk_id, layer_id, experts in stream_layers(path):
        if prev_chunk is None or chunk_id != prev_chunk:
            prev_chunk, prev_layer, prev_exps = chunk_id, layer_id, experts
            continue
        if prev_exps is not None and (layer_id == prev_layer + 1):
            key = (prev_layer, layer_id)
            for i in prev_exps:
                sources_by_pair[key].add(i)
                for j in experts:
                    targets_by_pair[key].add(j)
                    edges_by_pair[key][(i, j)] += 1
        prev_chunk, prev_layer, prev_exps = chunk_id, layer_id, experts

    return edges_by_pair, sources_by_pair, targets_by_pair

# ----------------------- Fast accumulator -----------------------

def accumulate_layerpair_edges_fast(path: str, max_expert_id: int):
    E = max_expert_id + 1
    mats_by_pair: Dict[Tuple[int,int], np.ndarray] = {}
    rows_by_pair: Dict[Tuple[int,int], Set[int]] = {}
    cols_by_pair: Dict[Tuple[int,int], Set[int]] = {}

    def get_mat(pair: Tuple[int,int]) -> np.ndarray:
        if pair not in mats_by_pair:
            mats_by_pair[pair] = np.zeros((E, E), dtype=np.uint32)
            rows_by_pair[pair] = set()
            cols_by_pair[pair] = set()
        return mats_by_pair[pair]

    prev_layer: Optional[int] = None
    prev_exps: Optional[List[int]] = None

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s[0] == '%' and all(c == '%' for c in s):
                prev_layer, prev_exps = None, None
                continue

            i = s.find(':')
            if i < 0:
                continue
            try:
                layer_id = int(s[:i].strip())
            except ValueError:
                continue

            exps: List[int] = []
            for tok in s[i+1:].split(','):
                tok = tok.strip()
                if not tok:
                    continue
                try:
                    v = int(tok)
                except ValueError:
                    continue
                if 0 <= v <= max_expert_id:
                    exps.append(v)

            if prev_layer is None:
                prev_layer, prev_exps = layer_id, exps
                continue

            if layer_id == prev_layer + 1 and prev_exps:
                key = (prev_layer, layer_id)
                M = get_mat(key)

                a = np.fromiter(prev_exps, dtype=np.int64)
                b = np.fromiter(exps, dtype=np.int64)
                rows_by_pair[key].update(prev_exps)
                cols_by_pair[key].update(exps)
                np.add.at(M, (a[:, None], b[None, :]), 1)

            prev_layer, prev_exps = layer_id, exps

    return mats_by_pair, rows_by_pair, cols_by_pair

# ----------------------- Clustering orders -----------------------

def order_rows_cols(M: np.ndarray, cluster_rows: bool, cluster_cols: bool, method: str = "fiedler") -> Tuple[List[int], List[int]]:
    r, c = M.shape
    row_order = list(range(r))
    col_order = list(range(c))

    def _safe_normalize_rows(X: np.ndarray) -> np.ndarray:
        s = X.sum(axis=1, keepdims=True)
        return np.divide(X, s, out=np.zeros_like(X, dtype=float), where=s != 0)

    def _cosine_row_sim(X: np.ndarray) -> np.ndarray:
        Xn = X.astype(float, copy=True)
        norms = np.linalg.norm(Xn, axis=1, keepdims=True)
        np.divide(Xn, norms, out=Xn, where=norms != 0)
        return Xn @ Xn.T

    def _fiedler_order(S: np.ndarray) -> List[int]:
        S = np.maximum(S, 0.0)
        S = (S + S.T) * 0.5
        d = S.sum(axis=1)
        L = np.diag(d) - S
        try:
            w, V = np.linalg.eigh(L)
        except np.linalg.LinAlgError:
            return list(range(S.shape[0]))
        if len(w) < 2:
            return list(range(S.shape[0]))
        idx = np.argsort(w)
        f = V[:, idx[1]] if len(idx) > 1 else V[:, idx[0]]
        return list(np.argsort(f))

    def _pca1_order(X: np.ndarray) -> List[int]:
        Xc = X - X.mean(axis=0, keepdims=True)
        try:
            U, S, _ = np.linalg.svd(Xc, full_matrices=False)
            pc1_scores = U[:, 0] * S[0] if S.size > 0 else np.zeros(X.shape[0])
        except np.linalg.LinAlgError:
            pc1_scores = np.zeros(X.shape[0])
        return list(np.argsort(pc1_scores))

    P = _safe_normalize_rows(M)

    if cluster_rows and r > 1:
        row_order = _fiedler_order(_cosine_row_sim(P)) if method == "fiedler" else _pca1_order(P)

    if cluster_cols and c > 1:
        PT = _safe_normalize_rows(M.T)
        col_order = _fiedler_order(_cosine_row_sim(PT)) if method == "fiedler" else _pca1_order(PT)

    return row_order, col_order

# ----------------------- Plotting -----------------------

def _downsample_step(n_labels: int, max_ticks: int) -> int:
    return max(1, n_labels // max(1, max_ticks))

def plot_matrix(out_dir: Path,
                M: np.ndarray,
                row_labels: List[int],
                col_labels: List[int],
                title: str,
                fname: str,
                max_ticks: int = 20,
                fontsize: int = 7) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.figure()
    im = plt.imshow(M, aspect="auto")
    plt.colorbar(im)

    sx = _downsample_step(len(col_labels), max_ticks)
    sy = _downsample_step(len(row_labels), max_ticks)

    plt.xticks(range(0, len(col_labels), sx), [str(x) for x in col_labels[::sx]], rotation=90, fontsize=fontsize)
    plt.yticks(range(0, len(row_labels), sy), [str(y) for y in row_labels[::sy]], fontsize=fontsize)

    plt.title(title)
    plt.tight_layout()
    out_path = out_dir / fname
    plt.savefig(out_path)
    plt.close()
    return out_path

# ----------------------- Build matrices -----------------------

def build_matrices_dict(path: str):
    edges_by_pair, sources_by_pair, targets_by_pair = accumulate_layerpair_edges_dict(path)
    mats: Dict[Tuple[int,int], Tuple[np.ndarray, List[int], List[int]]] = {}
    for (t, u), edges in edges_by_pair.items():
        rows = sorted(sources_by_pair[(t, u)])
        cols = sorted(targets_by_pair[(t, u)])
        if not rows or not cols:
            continue
        r_idx = {n: i for i, n in enumerate(rows)}
        c_idx = {n: j for j, n in enumerate(cols)}
        M = np.zeros((len(rows), len(cols)), dtype=float)
        for (i, j), w in edges.items():
            ri = r_idx.get(i)
            cj = c_idx.get(j)
            if ri is not None and cj is not None:
                M[ri, cj] += float(w)
        mats[(t, u)] = (M, rows, cols)
    return mats

def build_matrices_fast(path: str, max_id: int):
    mats_by_pair, rows_by_pair, cols_by_pair = accumulate_layerpair_edges_fast(path, max_expert_id=max_id)
    mats: Dict[Tuple[int,int], Tuple[np.ndarray, List[int], List[int]]] = {}
    for (t, u), Mfull in mats_by_pair.items():
        rows = sorted(rows_by_pair[(t, u)])
        cols = sorted(cols_by_pair[(t, u)])
        if not rows or not cols:
            continue
        Ms = Mfull[np.array(rows)[:, None], np.array(cols)[None, :]].astype(float, copy=False)
        mats[(t, u)] = (Ms, rows, cols)
    return mats

# ----------------------- Orchestration -----------------------

def parse_pairs_list(pairs: str) -> Optional[Set[Tuple[int,int]]]:
    pairs = (pairs or "").strip()
    if not pairs:
        return None
    out: Set[Tuple[int,int]] = set()
    for item in pairs.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            t_str, u_str = item.split(":")
            out.add((int(t_str), int(u_str)))
        except Exception:
            continue
    return out or None

def expected_outputs(mats_keys: List[Tuple[int,int]], out_dir: Path, probs: bool, only_probs: bool, order: str, cluster: str) -> List[Path]:
    outs: List[Path] = []
    for (t,u) in mats_keys:
        if probs or only_probs:
            outs.append(out_dir / f"clustered_L{t}_to_L{u}_probs_{order}_{cluster}.png")
            if probs and not only_probs:
                outs.append(out_dir / f"clustered_L{t}_to_L{u}_counts_like_{order}_{cluster}.png")
        else:
            outs.append(out_dir / f"clustered_L{t}_to_L{u}_counts_{order}_{cluster}.png")
    return outs

def run(cfg: dict) -> int:
    path = str(cfg["input_path"])
    input_path = os.path.expandvars(path)
    if "$" in input_path:
        raise RuntimeError(f"Unresolved env var in input_path: {path}")

    if not os.path.exists(input_path):
        raise FileNotFoundError(input_path)
    out_dir = Path(str(cfg.get("out_dir", "outputs")))
    method = str(cfg.get("method", "dict"))
    max_id = int(cfg.get("max_id", 63))
    probs = bool(cfg.get("probs", False))
    only_probs = bool(cfg.get("only_probs", False))
    cluster = str(cfg.get("cluster", "both"))
    order = str(cfg.get("order", "fiedler"))
    pairs = str(cfg.get("pairs", ""))
    max_ticks = int(cfg.get("max_ticks", 20))
    fontsize = int(cfg.get("fontsize", 7))
    force = bool(cfg.get("force", False))

    only_pairs = parse_pairs_list(pairs)

    if method == "dict":
        mats = build_matrices_dict(input_path)
    elif method == "fast":
        mats = build_matrices_fast(input_path, max_id=max_id)
    else:
        raise ValueError("expert_statistics.method must be 'dict' or 'fast'")

    if not mats:
        print("[expert_statistics] No consecutive layer pairs found.")
        return 0

    # filter to only_pairs (if provided) BEFORE skip check
    if only_pairs:
        mats = {k:v for k,v in mats.items() if k in only_pairs}
        if not mats:
            print("[expert_statistics] No matrices matched requested --pairs.")
            return 0

    keys = sorted(mats.keys())

    # Skip whole step if outputs already exist
    if not force:
        outs = expected_outputs(keys, out_dir, probs=probs, only_probs=only_probs, order=order, cluster=cluster)
        if outs and all(p.exists() for p in outs):
            print("[expert_statistics] Outputs already exist; skipping.")
            return 0

    for (t, u) in keys:
        M, rows, cols = mats[(t, u)]
        M_counts = M
        row_sums = M_counts.sum(axis=1, keepdims=True)
        M_probs = np.divide(M_counts, row_sums, out=np.zeros_like(M_counts), where=row_sums != 0)

        base = M_probs if (probs or only_probs) else M_counts

        cluster_rows = cluster in ("rows", "both")
        cluster_cols = cluster in ("cols", "both")

        r_ord, c_ord = order_rows_cols(base, cluster_rows, cluster_cols, method=order)
        Mr = base[np.array(r_ord)[:, None], np.array(c_ord)[None, :]]
        row_labs = [rows[i] for i in r_ord]
        col_labs = [cols[j] for j in c_ord]

        kind = "probabilities" if (probs or only_probs) else "counts"
        clabel = {"none":"no clustering","rows":"row clustering","cols":"column clustering","both":"row+column clustering"}[cluster]
        title = f"Cumulative {kind} (Layer {t} → {u})\n{clabel} via {order}"

        suffix = f"{'probs' if (probs or only_probs) else 'counts'}_{order}_{cluster}"
        fname = f"clustered_L{t}_to_L{u}_{suffix}.png"
        outp = plot_matrix(out_dir, Mr, row_labs, col_labs, title, fname, max_ticks=max_ticks, fontsize=fontsize)
        print(f"[expert_statistics] Saved: {outp}")

        if probs and not only_probs:
            Mr_counts = M_counts[np.array(r_ord)[:, None], np.array(c_ord)[None, :]]
            title2 = f"Cumulative counts (Layer {t} → {u})\nordered like {order} on probabilities"
            fname2 = f"clustered_L{t}_to_L{u}_counts_like_{order}_{cluster}.png"
            outp2 = plot_matrix(out_dir, Mr_counts, row_labs, col_labs, title2, fname2, max_ticks=max_ticks, fontsize=fontsize)
            print(f"[expert_statistics] Saved: {outp2}")

    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="pipeline.yaml")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        y = yaml.safe_load(f) or {}

    cfg = y.get("expert_statistics")
    if not isinstance(cfg, dict):
        raise SystemExit("YAML missing expert_statistics section")

    return run(cfg)

if __name__ == "__main__":
    raise SystemExit(main())
