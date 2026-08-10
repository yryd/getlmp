#!/usr/bin/env python
"""导出层·ReaxFF：坐标/元素 → 极简 data.lmp（Masses + Atoms，无键项）。

ReaxFF 的 LAMMPS 约定：
- atom_style charge（电荷列存在，QEq 模拟中计算，初始 0）
- 原子类型号 = reax_elements 中的元素顺序（须与 ffield.reax 的元素顺序一致）
- 无 Bonds/Angles/Dihedrals/Impropers 段（ReaxFF 键级由 pair_style reax/c 计算）
"""
from __future__ import annotations

import os

from check_lammps_data import parse_data

# 元素质量（g/mol，IUPAC 常规值；RDKit 不可用时兜底）
_ELEM_MASS = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38, 'Br': 79.904,
    'I': 126.90,
}


def build_data_reaxff(cfg, workdir: str, blocks: list, atom_info: list[dict],
                      box: list | None = None) -> dict:
    """ReaxFF 极简 data.lmp。

    blocks: list[type_index] -> list[分子拷贝] -> list[(elem, x, y, z)]
    atom_info: 每分子类型 {'elements': [...], 'natom': N}
    box: LAMMPS 顺序 [xlo,xhi,ylo,yhi,zlo,zhi]（None 时由坐标推算）
    """
    out_lmp = os.path.join(workdir, cfg.output)
    elem2type = {e: i + 1 for i, e in enumerate(cfg.reax_elements)}

    # 统计各类型用到的元素（Masses 段只列实际出现的）
    used_elems: list[str] = []
    for info in atom_info:
        for e in info['elements']:
            if e not in used_elems:
                used_elems.append(e)
    for e in used_elems:
        if e not in elem2type:
            raise ValueError(f'ReaxFF: 元素 {e} 不在 reax_elements {cfg.reax_elements} 中；'
                             f'请按 ffield.reax 的元素顺序在 yaml 配置 reax_elements')

    # 组装原子行
    # atom_style charge（默认）→ 6 列: id type q x y z（无 molecule 属性，ATOMIC 型）
    # atom_style full          → 7 列: id mol-id type q x y z（带分子 ID 便于分组输出）
    atoms_lines = []
    natom = 0
    mol_id = 0
    full = cfg.reax_atom_style == 'full'
    for ti, block in enumerate(blocks):
        for mol in block:
            mol_id += 1
            for (elem, x, y, z) in mol:
                natom += 1
                t = elem2type[elem]
                if full:
                    atoms_lines.append(
                        f'{natom:8d} {mol_id:8d} {t:4d} {0.0:12.6f} '
                        f'{x:16.6f} {y:16.6f} {z:16.6f}')
                else:
                    atoms_lines.append(
                        f'{natom:8d} {t:4d} {0.0:12.6f} '
                        f'{x:16.6f} {y:16.6f} {z:16.6f}')

    if box is None:
        box = _auto_box(atoms_lines, full)
    xlo, xhi, ylo, yhi, zlo, zhi = box
    if xhi <= xlo or yhi <= ylo or zhi <= zlo:
        raise ValueError(f'ReaxFF: 非法 box {box}')

    n_types = len(used_elems)
    lines = [
        f'LAMMPS data file from getlmp (ReaxFF, {cfg.name})',
        '',
        f'{natom:8d} atoms',
        f'{n_types:8d} atom types',
        '',
        f'{xlo:16.6f} {xhi:16.6f} xlo xhi',
        f'{ylo:16.6f} {yhi:16.6f} ylo yhi',
        f'{zlo:16.6f} {zhi:16.6f} zlo zhi',
        '',
        'Masses',
        '',
    ]
    for e in sorted(used_elems, key=lambda x: elem2type[x]):
        lines.append(f'{elem2type[e]:4d} {_mass(e):12.4f}  # {e}')
    lines += ['', 'Atoms  # charge', '']
    lines += atoms_lines
    lines += ['']
    with open(out_lmp, 'w') as f:
        f.write('\n'.join(lines))

    info = {
        'natom': natom, 'nbond': 0, 'nangle': 0, 'ndihedral': 0, 'nimproper': 0,
        'total_charge': 0.0, 'box': [xlo, xhi, ylo, yhi, zlo, zhi],
        'atom_types': n_types,
        'data_lmp': out_lmp,
    }
    check = parse_data(out_lmp)
    print(f'  data.lmp (ReaxFF): {info["natom"]} atoms / {n_types} atom types / '
          f'无键项 / charge=0.0 (QEq 模拟中算) / box={box}')
    return {'info': info, 'check': check, 'data_lmp': out_lmp}


def _mass(elem: str) -> float:
    try:
        from rdkit.Chem import PeriodicTable
        return float(PeriodicTable.GetAtomicWeight(elem))
    except Exception:
        return _ELEM_MASS.get(elem, 0.0)


def _auto_box(atoms_lines: list[str], full: bool) -> list:
    """从坐标推算盒（各方向 max-min + 2 Å padding）。"""
    xs, ys, zs = [], [], []
    # charge: id type q x y z → 坐标列 3,4,5；full: id mol type q x y z → 坐标列 4,5,6
    xc, yc, zc = (4, 5, 6) if full else (3, 4, 5)
    for ln in atoms_lines:
        p = ln.split()
        xs.append(float(p[xc])); ys.append(float(p[yc])); zs.append(float(p[zc]))
    pad = 2.0
    return [min(xs) - pad, max(xs) + pad,
            min(ys) - pad, max(ys) + pad,
            min(zs) - pad, max(zs) + pad]
