#!/usr/bin/env python
"""导出层核心：ParmEd Structure(prmtop) → LAMMPS data.lmp (atom_style full)

设计要点：
- 用 parmed 读 prmtop/inpcrd，获取全部力场参数（parmed 4.x 已移除 lammpsdata，
  故自写 LAMMPS data 输出；格式参考 MSTK lmpexporter + LAMMPS 规范）
- 拓扑段原子索引 ×3 编码由 parmed 内部处理；improper 以 dihedral.improper=True
  标记存在 s.dihedrals 中（s.impropers 为空），需手动筛出
- 二面角：Amber per/phase 项 → LAMMPS dihedral_style fourier 的 K1-K4
- improper：Amber per/phase → LAMMPS improper_style cvff 的 K, d, n
- 单位：real（Å, kcal/mol, e, g/mol），与 Amber 原生一致，无需换算
"""
from __future__ import annotations

from collections import defaultdict
from typing import List

import numpy as np
import parmed as pmd

RAD2DEG = 180.0 / np.pi


# ---------------------------------------------------------------- 工具函数

def _canon_dihedral_key(t1, t2, t3, t4):
    """二面角类型规范化：a-b-c-d 与 d-c-b-a 视为同一类型（LAMMPS 方向不敏感）"""
    fwd = (t1, t2, t3, t4)
    rev = (t4, t3, t2, t1)
    return fwd if fwd <= rev else rev


def _fourier_item(per, phase_deg, pk):
    """单个 Amber per/phase 项 → LAMMPS fourier 项 (K, n, d)。

    Amber:           E = pk·(1 + cos(nφ − γ)), γ ∈ {0°, 180°}
    LAMMPS fourier:  E = K·(1 + cos(nφ − d))   （逐项同构，直接映射）
    返回 (K, n, d)：K 能量(kcal/mol)，n 整数周期，d 度数。
    """
    n = int(round(per))
    if n < 0:
        raise ValueError(f'不支持 per={per}（LAMMPS fourier 要求 n≥0 整数）')
    phase = abs(float(phase_deg)) % 360.0
    if abs(phase) < 0.5 or abs(phase - 360.0) < 0.5:
        d = 0.0
    elif abs(phase - 180.0) < 0.5:
        d = 180.0
    else:
        raise ValueError(f'不支持 phase={phase_deg}（仅 0/180 可转 fourier）')
    return pk, n, d


def _cvff_improper(per, phase_deg, pk):
    """Amber improper per/phase → cvff 的 (K, d, n)。E = K(1 + d·cos(nφ))"""
    n = int(round(per))
    phase = abs(float(phase_deg)) % 360.0
    if abs(phase) < 0.5 or abs(phase - 360.0) < 0.5:
        d = 1.0
    elif abs(phase - 180.0) < 0.5:
        d = -1.0
    else:
        raise ValueError(f'不支持 phase={phase_deg}')
    return pk, d, n


# ---------------------------------------------------------------- 主转换

def prmtop_to_lammps(prmtop: str, inpcrd: str, out_lmp: str,
                     buffer: float = 3.8, box: list | None = None) -> dict:
    """prmtop/inpcrd → data.lmp。返回统计信息 dict。

    buffer: 无盒体系自动生成盒的真空层厚度（Å）
    box:    显式指定盒 [xlo, xhi, ylo, yhi, zlo, zhi]（优先于 buffer）
    """
    s = pmd.load_file(prmtop, inpcrd)

    # ---- 1. 原子类型 & Masses / Pair Coeffs -------------------------------
    atom_type_names: List[str] = []   # 类型名 -> LAMMPS type id (1-based)
    type_id = {}
    for a in s.atoms:
        # 4-site 水模型（OPC/TIP4P*）的虚拟位点（EP）：LAMMPS 无对应位置，
        # 需要专门导出链路（二期）；本期直接报错避免静默产错 data。
        if isinstance(a, pmd.ExtraPoint) or a.type == 'EP':
            raise RuntimeError(
                f'体系含 4-site 水虚拟位点（原子 {a.name!r}, type={a.type!r}）：'
                f'LAMMPS data 导出暂不支持 EP 位点（二期 OPC/TIP4P* 预留）。'
                f'请改用 3-site 水模型（tip3p / spce / opc3）')
        if a.type not in type_id:
            type_id[a.type] = len(atom_type_names) + 1
            atom_type_names.append(a.type)

    # LJ 参数：LJ_depth[type_idx-1]=ε, LJ_radius[type_idx-1]=Rmin/2
    # σ = 2*R / 2^(1/6) = R * 2^(5/6)
    lj_params = {}   # type name -> (epsilon, sigma)
    for tname in atom_type_names:
        idx = s.LJ_types[tname]
        eps = s.LJ_depth[idx - 1]
        rmin_half = s.LJ_radius[idx - 1]
        sigma = 2.0 * rmin_half / (2.0 ** (1.0 / 6.0))
        lj_params[tname] = (eps, sigma)

    # ---- 2. 键 / 角 类型 ------------------------------------------------
    bond_type_map = {}   # (k, req) -> id
    bond_types = []
    for b in s.bonds:
        key = (round(b.type.k, 6), round(b.type.req, 6))
        if key not in bond_type_map:
            bond_type_map[key] = len(bond_types) + 1
            bond_types.append(key)

    angle_type_map = {}  # (k, theteq) -> id
    angle_types = []
    for ag in s.angles:
        key = (round(ag.type.k, 6), round(ag.type.theteq, 6))
        if key not in angle_type_map:
            angle_type_map[key] = len(angle_types) + 1
            angle_types.append(key)

    # ---- 3. 二面角类型（按规范化 4 元组分组，同类型按 (per,phase) 去重后合并）---
    dih_type_map = {}    # canon key -> id
    dih_fourier = []     # id-1 -> [(K, n, d), ...] 多项列表
    dih_seen = defaultdict(set)   # canon key -> {(per, phase)}
    for d in s.dihedrals:
        if d.improper:
            continue     # improper 单独处理
        t = d.type
        key = _canon_dihedral_key(d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
        per_phase = (round(t.per, 3), round(abs(float(t.phase)) % 360.0, 3))
        if per_phase in dih_seen[key]:
            continue     # 同一物理二面角的重复实例，跳过（参数一致）
        if key not in dih_type_map:
            dih_type_map[key] = len(dih_fourier) + 1
            dih_fourier.append([])
        dih_seen[key].add(per_phase)
        item = _fourier_item(t.per, t.phase, t.phi_k)
        dih_fourier[dih_type_map[key] - 1].append(item)

    # ---- 4. improper 类型 ------------------------------------------------
    imp_type_map = {}    # canon key -> id
    imp_cvff = []        # id-1 -> (K, d, n)
    for d in s.dihedrals:
        if not d.improper:
            continue
        t = d.type
        key = _canon_dihedral_key(d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
        if key not in imp_type_map:
            imp_type_map[key] = len(imp_cvff) + 1
            imp_cvff.append(_cvff_improper(t.per, t.phase, t.phi_k))

    # ---- 5. 盒 ----------------------------------------------------------
    if box is not None:
        xlo, xhi, ylo, yhi, zlo, zhi = box
    else:
        coords = np.array([[a.xx, a.xy, a.xz] for a in s.atoms])
        lo = coords.min(axis=0) - buffer
        hi = coords.max(axis=0) + buffer
        xlo, ylo, zlo = lo
        xhi, yhi, zhi = hi

    # ---- 6. 组装 data 文本 -----------------------------------------------
    L = []
    L.append('LAMMPS data file from getlmp (prmtop -> data.lmp)')
    L.append('')
    natom, nbond, nangle = len(s.atoms), len(s.bonds), len(s.angles)
    ndihedral = sum(1 for d in s.dihedrals if not d.improper)
    nimproper = sum(1 for d in s.dihedrals if d.improper)
    L.append(f'{natom} atoms')
    L.append(f'{nbond} bonds')
    L.append(f'{nangle} angles')
    L.append(f'{ndihedral} dihedrals')
    L.append(f'{nimproper} impropers')
    L.append('')
    L.append(f'{len(atom_type_names)} atom types')
    L.append(f'{len(bond_types)} bond types')
    L.append(f'{len(angle_types)} angle types')
    L.append(f'{len(dih_fourier)} dihedral types')
    L.append(f'{len(imp_cvff)} improper types')
    L.append('')
    L.append(f'{xlo:.6f} {xhi:.6f} xlo xhi')
    L.append(f'{ylo:.6f} {yhi:.6f} ylo yhi')
    L.append(f'{zlo:.6f} {zhi:.6f} zlo zhi')
    L.append('')
    L.append('Masses')
    L.append('')
    for tname in atom_type_names:
        mass = next(a.mass for a in s.atoms if a.type == tname)
        L.append(f'{type_id[tname]:4d} {mass:12.5f}  # {tname}')
    L.append('')
    L.append('Pair Coeffs  # lj/cut')
    L.append('')
    for tname in atom_type_names:
        eps, sig = lj_params[tname]
        L.append(f'{type_id[tname]:4d} {eps:12.6f} {sig:12.6f}  # {tname}')
    L.append('')
    L.append('Bond Coeffs  # harmonic')
    L.append('')
    for i, (k, req) in enumerate(bond_types):
        L.append(f'{i + 1:4d} {k:12.4f} {req:10.4f}')
    L.append('')
    L.append('Angle Coeffs  # harmonic')
    L.append('')
    for i, (k, teq) in enumerate(angle_types):
        L.append(f'{i + 1:4d} {k:12.4f} {teq:10.4f}')
    L.append('')
    L.append('Dihedral Coeffs  # fourier')
    L.append('')
    for i, terms in enumerate(dih_fourier):
        # LAMMPS fourier 格式：typeID m K1 n1 d1 [K2 n2 d2 ...]，m=项数
        parts = [f'{i + 1:4d}', f'{len(terms):2d}']
        for K, n, d in terms:
            parts.append(f'{K:12.6f} {n:3d} {d:8.3f}')
        L.append(' '.join(parts))
    L.append('')
    L.append('Improper Coeffs  # cvff')
    L.append('')
    for i, (k, d, n) in enumerate(imp_cvff):
        # LAMMPS cvff 格式：typeID K d n，d 必须整数 ±1
        L.append(f'{i + 1:4d} {k:12.6f} {int(d):4d} {n:4d}')
    L.append('')
    L.append('Atoms  # full')
    L.append('')
    # 分子 id：按残基分组编号（Amber prmtop 无显式分子概念，残基即分子单位；
    # 同残基对象共享引用，用 (chain, number) 作键稳定；无残基信息时兜底 1）
    res_to_mol: dict = {}
    mol_ids: list[int] = []
    for a in s.atoms:
        r = a.residue
        if r is None:
            mol_ids.append(1)
            continue
        key = (r.chain, r.number)
        if key not in res_to_mol:
            res_to_mol[key] = len(res_to_mol) + 1
        mol_ids.append(res_to_mol[key])
    for i, a in enumerate(s.atoms):
        L.append(f'{i + 1:6d} {mol_ids[i]:6d} {type_id[a.type]:4d} {a.charge:12.6f} '
                 f'{a.xx:12.5f} {a.xy:12.5f} {a.xz:12.5f}  # {a.name} {a.residue.name}')
    L.append('')
    if nbond > 0:
        L.append('Bonds')
        L.append('')
        for i, b in enumerate(s.bonds):
            key = (round(b.type.k, 6), round(b.type.req, 6))
            L.append(f'{i + 1:6d} {bond_type_map[key]:4d} {b.atom1.idx + 1:6d} {b.atom2.idx + 1:6d}')
        L.append('')
    if nangle > 0:
        L.append('Angles')
        L.append('')
        for i, ag in enumerate(s.angles):
            key = (round(ag.type.k, 6), round(ag.type.theteq, 6))
            L.append(f'{i + 1:6d} {angle_type_map[key]:4d} '
                     f'{ag.atom1.idx + 1:6d} {ag.atom2.idx + 1:6d} {ag.atom3.idx + 1:6d}')
        L.append('')
    if ndihedral > 0:
        L.append('Dihedrals')
        L.append('')
        di = 0
        for d in s.dihedrals:
            if d.improper:
                continue
            di += 1
            key = _canon_dihedral_key(d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
            L.append(f'{di:6d} {dih_type_map[key]:4d} '
                     f'{d.atom1.idx + 1:6d} {d.atom2.idx + 1:6d} {d.atom3.idx + 1:6d} {d.atom4.idx + 1:6d}')
        L.append('')
    if nimproper > 0:
        # 注意：header 声明 0 impropers 时不允许出现 Impropers 段（哪怕为空），
        # LAMMPS 会报 "Invalid data file section: Impropers"
        L.append('Impropers')
        L.append('')
        ii = 0
        for d in s.dihedrals:
            if not d.improper:
                continue
            ii += 1
            key = _canon_dihedral_key(d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)
            # LAMMPS cvff improper: 中心原子是第 1 个；Amber 中心是第 2 个(j)
            # Amber 记录 (i,j,k,l) 中心=j；cvff 期望中心在前 → 写为 (j,i,k,l)
            a1, a2, a3, a4 = d.atom1, d.atom2, d.atom3, d.atom4
            L.append(f'{ii:6d} {imp_type_map[key]:4d} '
                     f'{a2.idx + 1:6d} {a1.idx + 1:6d} {a3.idx + 1:6d} {a4.idx + 1:6d}')
        L.append('')

    with open(out_lmp, 'w') as f:
        f.write('\n'.join(L))

    # ---- 7. 统计信息 -----------------------------------------------------
    info = {
        'natom': natom, 'nbond': nbond, 'nangle': nangle,
        'ndihedral': ndihedral, 'nimproper': nimproper,
        'natom_type': len(atom_type_names), 'nbond_type': len(bond_types),
        'nangle_type': len(angle_types), 'ndihedral_type': len(dih_fourier),
        'nimproper_type': len(imp_cvff),
        'total_charge': round(sum(a.charge for a in s.atoms), 6),
        'box': [xlo, xhi, ylo, yhi, zlo, zhi],
    }
    return info


# ---------------------------------------------------------------- CLI

def main():
    import sys
    if len(sys.argv) < 4:
        print('用法: python prmtop_to_lammps.py <prmtop> <inpcrd> <out.lmp> [buffer]')
        sys.exit(1)
    prmtop, inpcrd, out = sys.argv[1], sys.argv[2], sys.argv[3]
    buffer = float(sys.argv[4]) if len(sys.argv) > 4 else 3.8
    info = prmtop_to_lammps(prmtop, inpcrd, out, buffer=buffer)
    print(f'OK: {prmtop} -> {out}')
    for k, v in info.items():
        print(f'  {k}: {v}')


if __name__ == '__main__':
    main()
