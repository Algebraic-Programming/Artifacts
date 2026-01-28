# MoE Hypergraph Pipeline and ILP Optimization

This repository implements an end-to-end pipeline for analyzing Mixture-of-Experts (MoE) routing data, constructing hypergraphs and pairwise interaction graphs, and solving an Integer Linear Program (ILP) for expert partitioning under memory and MLA constraints.

The pipeline is fully configuration-driven via YAML and designed to be modular, reproducible, and extensible.

---

## Overview

The pipeline consists of the following high-level stages:

1. **Hyperedge extraction** from routing traces  
2. **Pair extraction** from hyperedges  
3. **Hypergraph construction** with guaranteed node coverage  
4. **ILP optimization** for expert partitioning  
5. **Expert statistics and analysis** (optional)

Each stage is implemented as a standalone script and orchestrated by a lightweight `main.py` runner using a YAML configuration file.

---

---

## Pipeline Execution

The pipeline is executed via:

```bash
pip install -r requirements.txt
python main.py --config pipeline.yaml
```

## Script Descriptions
1. extract_hyperedges.py
    - Parses MoE routing traces (described in .env.example)

    - Produces per-layer hyperedge CSVs

    - Each hyperedge represents a set of experts activated together

    - Outputs: 
        ```bash 
        hyperedges_k*/layer_XX.csv
        ```

2. extract_pairs.py

    - Converts hyperedges into pairwise expert interactions

    - Aggregates frequencies across hyperedges

    - Outputs: 
        ```bash 
        pairs/layer_XX_pairs.csv
        ```

3. hypergraph.py

    Core hypergraph construction logic:

    - Builds a global hypergraph from pair CSVs

    - Guarantees 100% node coverage (every node appears in ≥1 edge)

    - Supports greedy cover + percentile filtering

    - This module is used internally by the ILP solver.

4. run_ilp.py

    - Builds and solves the MLA v2 ILP model

    - Uses COPT as the solver backend

    - Supports:

        - Memory constraints

        - MLA replication limits

        - Optional next-layer MLA coupling

        - Warm starts from previous solutions

    - Outputs: 
        ```bash 
        solutions/partitioning_*.txt
        ```

5. ilp_parameters.py

    Defines the instanceParameter class:

    - Number of partitions

    - Memory bounds

    - MLA constraints

    - Objective selection

    - Solver limits and tolerances

6. expert_statistics.py

    Post-processing and analysis step:

    - Computes expert-to-expert transitions across layers

    - Builds cumulative transition matrices

    - Supports sparse or dense accumulation

    - Produces clustered heatmaps (counts or probabilities)

    - Outputs: 
        ```bash 
        outputs/*.png
        ```

## Environment Variables

Some scripts accept paths via environment variables.

Example .env:
```bash 
INPUT_PATH=/path/to/routing_trace.txt
```
