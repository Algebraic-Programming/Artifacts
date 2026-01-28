#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import yaml


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def main() -> int:
    # default to pipeline.yaml if not provided
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("pipeline.yaml")
    if not cfg_path.exists():
        print(f"Config not found: {cfg_path}")
        print("Usage: python main.py [pipeline.yaml]")
        return 2

    cfg = load_yaml(cfg_path)

    steps = cfg.get("pipeline", {}).get("steps", [])
    if not steps:
        print("pipeline.steps missing/empty in yaml")
        return 1

    for step in steps:
        name = step.get("name", "unnamed")
        script = step.get("script")
        if not script:
            print(f"Step {name} missing script")
            return 1

        cmd = [sys.executable, script, "--config", str(cfg_path)]

        print(f"\n==> {name}")
        print("    " + " ".join(cmd))

        subprocess.run(cmd, check=True)

    print("\n✅ done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
