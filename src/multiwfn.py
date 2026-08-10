#!/usr/bin/env python
"""Multiwfn 封装层：RESP / RESP2 电荷拟合 + ESP 可视化产物。

交互序列（2026-08-11 实测 Multiwfn 2026.7.15 noGUI）：
- RESP 电荷：  主功能 7 → 18 → 1 → y → 0 → 0 → q      → {base}.chg
- vtx.pdb：    主功能 12 → 0 → 6                       → vtx.pdb（闭合密度 0.001
                等值面顶点，B 因子字段存 ESP，单位 kcal/mol，VMD Beta 着色）
- 严格 ESP cube：主功能 12 → 0 → 13                    → mapfunc.cub（QM 电子云积分）
- RESP2：两次 RESP（气相+溶剂）后按 δ 混合（calcRESP.sh 逻辑）→ RESP2.chg

Multiwfn 交互式输入用 stdin 管道喂入；退出时可能报 Fortran I/O 错误
（exit 59，stdin 耗尽），只要目标产物文件已生成即视为成功。
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys

# 默认候选路径（config 未指定时按序探测）
DEFAULT_CANDS = [
    os.path.expanduser('~/packages/soft/multiwfn/Multiwfn_noGUI'),
    os.path.expanduser('~/packages/soft/multiwfn/*/Multiwfn_noGUI'),
    '/home/yryd/packages/soft/multiwfn/Multiwfn_noGUI',
    '/home/yryd/packages/soft/multiwfn/*/Multiwfn_noGUI',
]

# Multiwfn RESP 标准输入（官方 calcRESP.sh）
_RESP_INPUT = '7\n18\n1\ny\n0\n0\nq\n'
# 表面分析：12→0 开始分析 → 6 导出 vtx.pdb → -1 返回
_ESP_SURFACE_INPUT = '12\n0\n6\n-1\n0\nq\n'
# 表面分析：12→0 开始分析 → 13 导出 mapfunc.cub → -1 返回
_ESP_CUBE_INPUT = '12\n0\n13\n-1\n0\nq\n'


def find_multiwfn(cfg_path: str = '') -> str:
    """探测 Multiwfn_noGUI 可执行路径。

    顺序：config 显式路径 → PATH → 默认候选（含 glob 通配）。
    找不到时报错并给出提示。
    """
    if cfg_path:
        if os.path.isfile(cfg_path):
            return cfg_path
        raise FileNotFoundError(f'multiwfn_path 不存在: {cfg_path}')
    # PATH
    for name in ('Multiwfn_noGUI', 'multiwfn'):
        p = shutil.which(name)
        if p:
            return p
    # 默认候选（含 glob）
    for cand in DEFAULT_CANDS:
        hits = glob.glob(cand)
        if hits and os.path.isfile(hits[0]):
            return hits[0]
    raise FileNotFoundError(
        '找不到 Multiwfn_noGUI。请在 yaml 配置 qm.multiwfn_path 指定路径，'
        '或安装到 ~/packages/soft/multiwfn/ 下')


def _run_multiwfn(mwfn: str, wfn_file: str, stdin_text: str,
                  cwd: str, timeout: int = 900) -> str:
    """以 stdin 管道喂输入运行 Multiwfn，返回 stdout 文本。

    退出码非 0（如 59 = stdin 耗尽）不视为失败——由调用方检查产物文件。
    """
    r = subprocess.run(
        [mwfn, os.path.abspath(wfn_file), '-ispecial', '1'],
        cwd=cwd, input=stdin_text, capture_output=True, text=True,
        timeout=timeout,
    )
    if r.returncode != 0 and r.returncode != 59:
        print(f'  [warn] Multiwfn 退出码 {r.returncode}（可能仍已产出文件）',
              file=sys.stderr)
    return r.stdout


def _parse_chg(path: str) -> list[float]:
    """解析 Multiwfn .chg：每行 name X Y Z charge → [charge,...]。"""
    chgs = []
    with open(path) as f:
        for ln in f:
            parts = ln.split()
            if len(parts) >= 5:
                try:
                    chgs.append(float(parts[4]))
                except ValueError:
                    continue
    if not chgs:
        raise RuntimeError(f'Multiwfn .chg 未解析到电荷: {path}')
    return chgs


def resp_from_fch(fch: str, workdir: str = '.',
                  multiwfn_path: str = '') -> tuple[list[float], str]:
    """对含 MK ESP 数据的 .fch 拟合 RESP 电荷。

    返回 (charges, chg_path)。chg_path 为 Multiwfn 输出的 {base}.chg。
    """
    mwfn = find_multiwfn(multiwfn_path)
    base = os.path.splitext(os.path.basename(fch))[0]
    chg_path = os.path.join(workdir, base + '.chg')
    if os.path.exists(chg_path):
        os.remove(chg_path)
    _run_multiwfn(mwfn, fch, _RESP_INPUT, cwd=workdir)
    if not os.path.exists(chg_path):
        raise RuntimeError(f'Multiwfn RESP 未生成电荷文件: {chg_path}')
    return _parse_chg(chg_path), chg_path


def resp2_from_fch(fch_gas: str, fch_solv: str, delta: float = 0.5,
                   workdir: str = '.', multiwfn_path: str = '') -> list[float]:
    """RESP2：气相与溶剂各拟合 RESP，再按 q = (1-δ)*q_gas + δ*q_solv 混合。

    对应官方 calcRESP.sh 逻辑（Multiwfn 无单命令 RESP2，两次 RESP + awk）。
    """
    q_gas, _ = resp_from_fch(fch_gas, workdir, multiwfn_path)
    q_solv, _ = resp_from_fch(fch_solv, workdir, multiwfn_path)
    if len(q_gas) != len(q_solv):
        raise RuntimeError(f'RESP2 原子数不一致: gas={len(q_gas)} solv={len(q_solv)}')
    return [(1 - delta) * a + delta * b for a, b in zip(q_gas, q_solv)]


def esp_surface_pdb(fch: str, out_pdb: str, workdir: str = '.',
                    multiwfn_path: str = '') -> str:
    """生成闭合 ESP 等值面顶点文件 vtx.pdb（B 因子 = ESP kcal/mol）。

    密度 0.001 闭合等值面（Multiwfn 主功能 12 默认）+ 表面 ESP 映射。
    返回输出路径。
    """
    mwfn = find_multiwfn(multiwfn_path)
    _run_multiwfn(mwfn, fch, _ESP_SURFACE_INPUT, cwd=workdir)
    src = os.path.join(workdir, 'vtx.pdb')
    if not os.path.exists(src):
        raise RuntimeError(f'Multiwfn 未生成 vtx.pdb（ESP 曲面失败）: {src}')
    shutil.copy(src, out_pdb)
    return out_pdb


def esp_cube(fch: str, out_cub: str, workdir: str = '.',
             multiwfn_path: str = '') -> str:
    """生成严格 QM ESP cube（mapfunc.cub，QM 电子云积分，非点电荷近似）。

    返回输出路径。
    """
    mwfn = find_multiwfn(multiwfn_path)
    _run_multiwfn(mwfn, fch, _ESP_CUBE_INPUT, cwd=workdir)
    src = os.path.join(workdir, 'mapfunc.cub')
    if not os.path.exists(src):
        raise RuntimeError(f'Multiwfn 未生成 mapfunc.cub（ESP cube 失败）: {src}')
    shutil.copy(src, out_cub)
    return out_cub
