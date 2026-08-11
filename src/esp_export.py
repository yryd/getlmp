# -*- coding: utf-8 -*-
"""esp_export.py — 从 LAMMPS data.lmp 生成静电势 .cub（点电荷库仑近似）。

复刻 Multiwfn outcube 的 Gaussian cube 写出格式；ESP 值 = Σ q_i / |r - r_i|
（a.u.），即 Multiwfn cubesp 核项 nucesp 的模拟（核电荷 Z → 力场电荷 q）。
numpy 向量化网格计算，单分子毫秒级。

用法（模块内）：
    from esp_export import export_esp_cube
    export_esp_cube('data.lmp', 'out.cub', spacing=0.3, buffer=1.5)
"""
from __future__ import annotations

import os

import numpy as np

MASS2ELEM = {1.008: 'H', 4.0026: 'He', 6.94: 'Li', 9.0122: 'Be', 10.81: 'B',
             12.011: 'C', 14.007: 'N', 15.999: 'O', 18.998: 'F', 20.180: 'Ne',
             22.990: 'Na', 24.305: 'Mg', 26.982: 'Al', 28.085: 'Si', 30.974: 'P',
             32.06: 'S', 35.45: 'Cl', 39.948: 'Ar', 39.098: 'K', 40.078: 'Ca',
             55.845: 'Fe', 63.546: 'Cu', 65.38: 'Zn', 79.904: 'Br', 126.904: 'I'}
ELEM2NUM = {e: n for n, e in enumerate(
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
     'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'], 1)}


def parse_data_lmp(path: str) -> tuple[list, dict]:
    """解析 data.lmp → (atoms, masses)。
    atoms: [(x,y,z,q,type)]; masses: {type: mass}。"""
    atoms, masses = [], {}
    section = None
    with open(path) as f:
        for ln in f:
            ln = ln.split('#')[0].strip()
            if not ln or ln.isdigit():
                continue
            up = ln.upper()
            if up.startswith('MASSES'):
                section = 'masses'
                continue
            if up.startswith('ATOMS'):
                section = 'atoms'
                continue
            if up.startswith(('BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS',
                              'PAIR', 'BOND', 'ANGLE', 'DIHEDRAL', 'IMPROPER',
                              'VELOCITIES', 'UNITS', 'ATOM', 'HEADER', 'BY')):
                section = None
                continue
            parts = ln.split()
            if section == 'masses' and len(parts) >= 2:
                masses[int(parts[0])] = float(parts[1])
            elif section == 'atoms' and len(parts) >= 6:
                x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
                q = float(parts[-4])
                typ = int(parts[2]) if len(parts) >= 7 else int(parts[1])
                atoms.append((x, y, z, q, typ))
    if not atoms:
        raise RuntimeError(f'未解析到原子: {path}')
    return atoms, masses


def _elem_of_mass(m: float) -> str:
    best, bd = 'X', 1e9
    for mass, ele in MASS2ELEM.items():
        d = abs(m - mass)
        if d < bd:
            best, bd = ele, d
    return best if bd < 0.5 else 'X'


def compute_esp_cube(atoms: list, spacing: float = 0.3, buffer: float = 1.5,
                     min_r: float = 0.3):
    """numpy 向量化：在规则网格上算点电荷库仑 ESP（a.u.）。

    返回 (data, org, (nx,ny,nz))。data 形状 (nx,ny,nz)，单位 a.u.。
    """
    coords = np.array([[a[0], a[1], a[2]] for a in atoms], dtype=float)
    lo = coords.min(axis=0) - buffer
    hi = coords.max(axis=0) + buffer
    n = np.ceil((hi - lo) / spacing).astype(int) + 1
    org = lo

    # 网格坐标轴（1D）
    axes = [org[k] + np.arange(n[k]) * spacing for k in range(3)]
    # 3D 网格（广播，内存 ~nx*ny*nz*24B，单分子量级无压力）
    gx, gy, gz = np.meshgrid(*axes, indexing='ij')
    data = np.zeros_like(gx)
    for (ax, ay, az, q, _t) in atoms:
        r = np.sqrt((gx - ax) ** 2 + (gy - ay) ** 2 + (gz - az) ** 2)
        r = np.maximum(r, min_r)          # 原子核处防除零
        data += 0.52917721092 * q / r     # e/Å → a.u.
    return data, org, tuple(n.tolist())


def write_cube(path: str, atoms: list, masses: dict, data, org, n, spacing):
    """照 Multiwfn outcube 模板写 Gaussian cube。"""
    nx, ny, nz = n
    with open(path, 'w') as f:
        f.write(' ESP by getlmp (point-charge approx, a.u.)\n')
        f.write(f' Totally {nx*ny*nz} grid points\n')
        f.write(f'{len(atoms):5d}{org[0]:12.6f}{org[1]:12.6f}{org[2]:12.6f}\n')
        f.write(f'{nx:5d}{spacing:12.6f}{0.0:12.6f}{0.0:12.6f}\n')
        f.write(f'{ny:5d}{0.0:12.6f}{spacing:12.6f}{0.0:12.6f}\n')
        f.write(f'{nz:5d}{0.0:12.6f}{0.0:12.6f}{spacing:12.6f}\n')
        for (ax, ay, az, q, typ) in atoms:
            m = masses.get(typ, 1.0)
            znum = ELEM2NUM.get(_elem_of_mass(m), 1)
            f.write(f'{znum:5d}{float(znum):12.6f}{ax:12.6f}{ay:12.6f}{az:12.6f}\n')
        # 数据：x→y→z 顺序，每行 6 个值（1PE14.5E3 风格）
        for i in range(nx):
            for j in range(ny):
                vals = data[i, j, :]
                for s in range(0, nz, 6):
                    chunk = vals[s:s + 6]
                    f.write(''.join(f'{v:14.5E}' for v in chunk) + '\n')


def export_esp_cube(data_lmp: str, out_cub: str,
                    spacing: float = 0.3, buffer: float = 1.5) -> str:
    """一键：data.lmp → ESP cube。返回输出路径。"""
    atoms, masses = parse_data_lmp(data_lmp)
    data, org, n = compute_esp_cube(atoms, spacing, buffer)
    write_cube(out_cub, atoms, masses, data, org, n, spacing)
    return out_cub


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='data.lmp → ESP cube（点电荷近似）')
    ap.add_argument('data_lmp')
    ap.add_argument('out_cub', nargs='?', default=None)
    ap.add_argument('--spacing', type=float, default=0.3)
    ap.add_argument('--buffer', type=float, default=1.5)
    args = ap.parse_args()
    out = args.out_cub or (os.path.splitext(args.data_lmp)[0] + '_esp.cub')
    export_esp_cube(args.data_lmp, out, args.spacing, args.buffer)
    print(f'写出 {out}')
