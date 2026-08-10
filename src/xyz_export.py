#!/usr/bin/env python
"""导出层：体系标准 xyz 导出（流水线默认主产物，workdir/system.xyz）。

标准 xyz 格式：第 1 行原子数；第 2 行注释；随后每行 `元素 x y z`（Å）。
**原子顺序与 data.lmp 完全一致**（可视化/后续工具可直接对应）：

- GAFF2（单/多分子）：prmtop + inpcrd → 元素（原子名）+ 坐标。
  tleap loadpdb 会按 prepi 模板重排残基内原子顺序，packed.xyz 原序
  与 data.lmp 不一致，因此必须从 prmtop/inpcrd（tleap 后的顺序）导出。
- ReaxFF 单分子：分子层坐标块（元素 + 坐标已在手）。
- ReaxFF 多分子：packed.xyz（无 tleap，packmol 顺序即 data.lmp 顺序）。
"""
from __future__ import annotations

import os
import shutil

import parmed as pmd


def export_system_xyz(report: dict, out_path: str) -> int:
    """把最终体系坐标导出为标准 xyz，返回原子数。"""
    cfg = report['config']
    workdir = report['workdir']
    mols = report['molecules']
    multi = cfg.packmol.enabled

    if cfg.forcefield == 'reaxff':
        if multi:
            # packed.xyz 顺序 = data.lmp（ReaxFF 无 tleap 重排）
            src = os.path.join(workdir, 'packed.xyz')
            if not os.path.exists(src):
                raise RuntimeError(f'体系 xyz 源缺失: {src}')
            with open(src) as f:
                n = int(f.readline().split()[0])
            shutil.copyfile(src, out_path)
            return n
        # ReaxFF 单分子：分子层坐标块
        m0 = mols[0]
        return _write_xyz(out_path, m0['elements'], m0['coords'],
                          f'from getlmp ({cfg.forcefield}) {cfg.name}')

    # GAFF2：统一从 prmtop/inpcrd（tleap 重排后顺序 = data.lmp 顺序）
    sysr = report['system']
    elements = _prmtop_elements(sysr['prmtop'])
    coords = _read_inpcrd_coords(sysr['inpcrd'], len(elements))
    return _write_xyz(out_path, elements, coords,
                      f'from getlmp ({cfg.forcefield}) {cfg.name}')


def _prmtop_elements(prmtop: str) -> list[str]:
    """prmtop 原子名 → 元素（顺序 = prmtop 原子顺序 = data.lmp 顺序）。"""
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    s = pmd.load_file(prmtop)
    elements: list[str] = []
    for a in s.atoms:
        sym = ''.join(ch for ch in a.name if not ch.isdigit())
        z = pt.GetAtomicNumber(sym)
        if z == 0:
            raise RuntimeError(f'无法从原子名推断元素: {a.name!r} ({prmtop})')
        elements.append(pt.GetElementSymbol(z))
    return elements


def _write_xyz(path: str, elements: list[str], coords: list[list[float]],
               comment: str) -> int:
    lines = [str(len(elements)), comment]
    for e, c in zip(elements, coords):
        lines.append(f'{e:>2s} {c[0]:12.5f} {c[1]:12.5f} {c[2]:12.5f}')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return len(elements)


def _read_inpcrd_coords(inpcrd: str, natom: int) -> list[list[float]]:
    """读 Amber inpcrd 坐标（第 1 行标题，第 2 行 natom，之后每行 6 个 12.7f 值）。

    只取前 natom*3 个数值（防行 2 出现盒信息等额外字段）。
    """
    with open(inpcrd) as f:
        lines = f.read().splitlines()
    vals: list[float] = []
    for ln in lines[2:]:
        for v in ln.split():
            try:
                vals.append(float(v))
            except ValueError:
                continue
            if len(vals) >= natom * 3:
                break
        if len(vals) >= natom * 3:
            break
    if len(vals) < natom * 3:
        raise RuntimeError(f'inpcrd 坐标数不足: {len(vals)} < {natom * 3} ({inpcrd})')
    return [[vals[i], vals[i + 1], vals[i + 2]] for i in range(0, natom * 3, 3)]
