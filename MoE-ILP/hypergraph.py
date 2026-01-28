import re
import glob
import csv
import math
from pathlib import Path
from collections import Counter, defaultdict
from typing import Iterable, List, Dict, Set, Optional

MAX_NODES_PER_LAYER_DEFAULT = 64
MAX_LAYERS_DEFAULT = 26
MEM_EXPERT_IN_GB = 0.08
MEM_MLA_IN_GB = 0.533

class Hypergraph:
    def __init__(
        self,
        edges: Iterable[Iterable[int]],
        node_mem_weights: List[float],
        edge_weights: List[float],
        node_names: List[str],
        original_id: List[int],
        id_map: Dict[int, int],
    ):
        self.edges = [set(e) for e in edges]
        self.edge_weights = list(edge_weights)
        self.original_id = list(original_id)
        self.id_map = dict(id_map)
        self.node_names = list(node_names)
        self.node_mem_weights = list(node_mem_weights)

        self.num_nodes = len(self.original_id)
        self.num_edges = len(self.edges)

        if not (len(self.node_mem_weights) == len(self.node_names) == len(self.original_id)):
            raise ValueError("Node attribute lengths mismatch")
        if len(self.edge_weights) != self.num_edges:
            raise ValueError("edge_weights length mismatch")

    @staticmethod
    def _extract_layer_index(path: str) -> int:
        name = Path(path).name
        m = re.search(r"(?:layer|L)[-_]?(\d+)", name, flags=re.IGNORECASE)
        if not m:
            m = re.search(r"(\d+)", name)
        if not m:
            raise ValueError(f"Cannot determine layer index from {name}")
        return int(m.group(1))

    def ensure_nodes_exist(
        self,
        gids: Iterable[int],
        *,
        max_nodes_per_layer: int = MAX_NODES_PER_LAYER_DEFAULT,
        default_mem: float = MEM_EXPERT_IN_GB,
    ) -> None:
        for gid in gids:
            if gid in self.id_map:
                continue
            cid = len(self.original_id)
            self.original_id.append(gid)
            self.id_map[gid] = cid
            layer, local = divmod(gid, max_nodes_per_layer)
            self.node_names.append(f"L{layer}_n{local}")
            self.node_mem_weights.append(default_mem)
        self.num_nodes = len(self.original_id)

    @classmethod
    def from_multilayer_csvs_pairs(
        cls,
        pattern: str,
        layer_threshold: int,
        per_layer_top_q: float = 0.98,
        max_nodes_per_layer: int = MAX_NODES_PER_LAYER_DEFAULT,
        max_layers: int = MAX_LAYERS_DEFAULT,
        *,
        strict_layer_files: bool = True,
        debug: bool = True,
        enforce_full_cover: bool = True,
        **_ignored,
    ) -> "Hypergraph":
        """
        Build from pair CSVs and ensure ALL nodes (0..max_nodes_per_layer-1) in each included layer
        are present AND covered by at least one selected edge.

        Fixes layer indexing:
        - Filenames use 1-based layers: layer_01_pairs.csv -> layer1 = 1
        - Global id: gid = (layer1 - 1) * max_nodes_per_layer + local_id
        - Per-layer percentile cuts are keyed by layer1 (1-based)
        """
        files = sorted(glob.glob(pattern))
        if not files:
            raise FileNotFoundError(pattern)

        if strict_layer_files:
            valid_re = re.compile(r"^layer_\d+_pairs\.csv$", flags=re.IGNORECASE)
            files = [p for p in files if valid_re.match(Path(p).name)]
            if not files:
                raise FileNotFoundError(
                    "No files passed strict filename check (expected like 'layer_01_pairs.csv')."
                )

        # --- Read all pairs (raw), aggregate, and collect per-layer weights ---
        pair_counter: Counter[tuple[int, int]] = Counter()
        pair_layer: Dict[tuple[int, int], int] = {}      # (g1,g2) -> layer1
        layer_weights: Dict[int, List[float]] = defaultdict(list)
        layers_seen: Set[int] = set()

        for fp in files:
            layer1 = cls._extract_layer_index(fp)   # expects 1-based from filename
            if not (1 <= layer1 <= max_layers):
                continue
            if layer1 > layer_threshold:
                continue

            layers_seen.add(layer1)

            with open(fp, newline="", encoding="utf-8") as f:
                rdr = csv.DictReader(f)
                for r in rdr:
                    try:
                        a = int(r["expert_a"])
                        b = int(r["expert_b"])
                        w = float(r["frequency"])
                    except Exception:
                        continue

                    if not (0 <= a < max_nodes_per_layer) or not (0 <= b < max_nodes_per_layer):
                        continue
                    if a == b:
                        continue

                    ga = (layer1 - 1) * max_nodes_per_layer + a
                    gb = (layer1 - 1) * max_nodes_per_layer + b
                    g1, g2 = sorted((ga, gb))

                    pair_counter[(g1, g2)] += w
                    pair_layer[(g1, g2)] = layer1
                    layer_weights[layer1].append(w)

        if not layers_seen:
            raise ValueError("No layers included (check layer_threshold / filenames / max_layers).")

        # --- Make ALL nodes exist for included layers ---
        all_nodes: Set[int] = set()
        for layer1 in layers_seen:
            base = (layer1 - 1) * max_nodes_per_layer
            all_nodes.update(range(base, base + max_nodes_per_layer))

        original_id = sorted(all_nodes)
        id_map = {g: i for i, g in enumerate(original_id)}
        node_mem = [MEM_EXPERT_IN_GB] * len(original_id)
        node_names = [f"L{g // max_nodes_per_layer}_n{g % max_nodes_per_layer}" for g in original_id]

        # If there are no pairs at all, we cannot cover nodes with edges.
        if not pair_counter:
            raise ValueError("No pairs found in input; cannot build edges / cover nodes.")

        # --- Per-layer percentile cut (1-based layer keys) ---
        def percentile(vals: List[float], q: float) -> float:
            if not vals:
                return 0.0
            vals = sorted(vals)
            k = max(1, min(len(vals), int(math.ceil(q * len(vals)))))
            return float(vals[k - 1])

        qcut = {L: percentile(ws, per_layer_top_q) for L, ws in layer_weights.items()}

        # --- Initial selection: keep edges passing per-layer cutoff ---
        selected_pairs: Dict[tuple[int, int], float] = {}
        for (g1, g2), w in pair_counter.items():
            L = pair_layer[(g1, g2)]  # 1-based
            if w >= qcut.get(L, 0.0):
                selected_pairs[(g1, g2)] = float(w)

        # --- Enforce FULL COVER: add best edges to cover any uncovered node ---
        if enforce_full_cover:
            # Build adjacency from ALL raw pairs (not just selected)
            best_incident: Dict[int, tuple[tuple[int, int], float]] = {}
            for (g1, g2), w in pair_counter.items():
                # best edge incident to g1
                cur = best_incident.get(g1)
                if cur is None or w > cur[1]:
                    best_incident[g1] = ((g1, g2), float(w))
                # best edge incident to g2
                cur = best_incident.get(g2)
                if cur is None or w > cur[1]:
                    best_incident[g2] = ((g1, g2), float(w))

            # Compute covered nodes under current selected edges
            covered: Set[int] = set()
            for (g1, g2) in selected_pairs.keys():
                covered.add(g1)
                covered.add(g2)

            missing = set(original_id) - covered
            if debug:
                print(f"[pairs] initially selected edges={len(selected_pairs)}")
                print(f"[pairs] initially covered={len(covered)} / {len(original_id)}; missing={len(missing)}")

            # Add best edge for each missing node (must exist in raw data)
            added = 0
            still_missing = []
            for gid in sorted(missing):
                if gid in covered:
                    continue
                best = best_incident.get(gid)
                if best is None:
                    still_missing.append(gid)
                    continue
                (g1, g2), w = best
                if (g1, g2) not in selected_pairs:
                    selected_pairs[(g1, g2)] = float(w)
                    added += 1
                covered.add(g1)
                covered.add(g2)

            # Re-check
            missing2 = set(original_id) - covered
            if missing2:
                # This can happen only if some nodes never appear in ANY pair row.
                # In that case, you cannot cover them with pair-edges.
                raise ValueError(
                    f"Cannot achieve 100% cover: {len(missing2)} nodes never appear in any pair. "
                    f"Examples: {sorted(list(missing2))[:25]}"
                )

            if debug:
                print(f"[pairs] added {added} edges for coverage repair; final edges={len(selected_pairs)}")
                print(f"[pairs] final covered={len(covered)} / {len(original_id)} (100%)")

        # --- Finalize edges in COMPACT ids ---
        edges_global = [(g1, g2) for (g1, g2) in selected_pairs.keys()]
        edge_weights = [selected_pairs[(g1, g2)] for (g1, g2) in selected_pairs.keys()]
        edges = [{id_map[g1], id_map[g2]} for (g1, g2) in edges_global]

        return cls(edges, node_mem, edge_weights, node_names, original_id, id_map)


    def add_transitions(
        self,
        folder: str | Path,
        per_layer_top_q: float,
        weight: float,
        layer_threshold: int,
        filename_glob: str = "groups_*_to_*.csv",
        *,
        sources_col: str = "sources",
        weight_col: str = "weight",
        list_delim: str = ";",
    ) -> None:
        folder = Path(folder)
        files = sorted(folder.glob(filename_glob))
        if not files:
            raise FileNotFoundError(f"No CSV files matching '{filename_glob}' in: {folder}")

        layer_weights: Dict[int, List[float]] = defaultdict(list)

        def layer_from_name(p: Path) -> int:
            m = re.search(r"(\d+)", p.name)
            if not m:
                raise ValueError(f"Cannot extract layer id from filename: {p.name}")
            return int(m.group(1))

        for fp in files:
            L = layer_from_name(fp)
            if L > layer_threshold:
                continue
            with open(fp, newline="", encoding="utf-8") as f:
                rdr = csv.DictReader(f)
                for r in rdr:
                    w = float(r[weight_col])
                    layer_weights[L].append(w)

        def percentile(vals: List[float], q: float) -> float:
            if not vals:
                return 0.0
            vals = sorted(vals)
            k = max(1, min(len(vals), int(math.ceil(q * len(vals)))))
            return float(vals[k - 1])

        qcut = {L: percentile(ws, per_layer_top_q) for L, ws in layer_weights.items()}

        for fp in files:
            L = layer_from_name(fp)
            if L > layer_threshold:
                continue
            cutoff = qcut.get(L, 0.0) * weight
            with open(fp, newline="", encoding="utf-8") as f:
                rdr = csv.DictReader(f)
                for r in rdr:
                    gids = [(L - 1) * MAX_NODES_PER_LAYER_DEFAULT + int(x) for x in str(r[sources_col]).split(list_delim) if x.strip() != ""]
                    if not gids:
                        continue
                    self.ensure_nodes_exist(gids, max_nodes_per_layer=MAX_NODES_PER_LAYER_DEFAULT)
                    w = float(r[weight_col])
                    if w * weight >= cutoff:
                        edge = {self.id_map[g] for g in gids}
                        if len(edge) > 1:
                            self.edges.append(edge)
                            self.edge_weights.append(w)

        self.num_edges = len(self.edges)

    def save(self, filepath: str) -> None:
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["experts", "frequency"])
            for e, wt in zip(self.edges, self.edge_weights):
                w.writerow([";".join(map(str, sorted(e))), wt])

    def compact_to_original(self, cid: int) -> int:
        return self.original_id[cid]

    def original_to_compact(self, gid: int) -> int:
        return self.id_map[gid]
    
    def edge_to_layer(
    self,
    edge_index: int,
    ) -> int:
        """
        Return the 1-based layer index for the given edge.
        Assumes all nodes in an edge are from the same layer.
        """
        if not (0 <= edge_index < self.num_edges):
            raise IndexError(f"edge_index {edge_index} out of range [0, {self.num_edges-1}]")

        edge = self.edges[edge_index]
        if not edge:
            raise ValueError(f"Edge {edge_index} is empty.")

        # Pick any node in the edge
        nid = next(iter(edge))

        # If edges use compact ids, map to original; if they already store original ids,
        # fall back to identity (works when id_map is identity).
        if 0 <= nid < len(self.original_id):
            orig_gid = self.original_id[nid]
        else:
            raise ValueError(
                f"Cannot resolve node id {nid} from edge {edge_index} to an original global id."
            )

        return (orig_gid // MAX_NODES_PER_LAYER_DEFAULT) + 1
