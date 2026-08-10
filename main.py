#!/usr/bin/env python
"""getlmp：SMILES + yaml → LAMMPS data.lmp 一键工具链。

用法：
    python main.py input.yaml
    python main.py input.yaml --buffer 5.0

依赖（RDKit / parmed / antechamber / packmol / tleap）在 conda 环境内运行。
"""
from __future__ import annotations

import argparse
import os
import sys

# 让 src/ 下的模块可直接导入（零安装运行）
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from config import load_config  # noqa: E402
from pipeline import run_pipeline  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="getlmp",
        description="SMILES + yaml → LAMMPS data.lmp 一键工具链"
                    "（GAFF2/AM1-BCC 等）",
    )
    ap.add_argument("input", help="input.yaml 配置文件")
    ap.add_argument("--buffer", type=float, default=None,
                    help="单分子盒 padding（Å，默认 3.8）")
    ap.add_argument("--version", action="version", version="getlmp 0.1.0")
    args = ap.parse_args(argv)

    cfg = load_config(args.input)
    if args.buffer is not None:
        cfg.buffer = args.buffer
    report = run_pipeline(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
