#!/usr/bin/env python
"""体系层：单分子 tleap / 多分子 packmol 坐标组装 + tleap。

- 单分子（阶段 1）：mol2 + frcmod → prmtop/inpcrd
- 多分子（阶段 2）：packed.xyz 坐标 + 各类型 mol2 → 合并 PDB → tleap loadpdb → prmtop/inpcrd
"""
from __future__ import annotations

import os
import subprocess
import sys

import parmed as pmd

from config import Config


def _run(cmd: list[str], cwd: str, desc: str) -> None:
    print(f'  [run] {" ".join(cmd)}', file=sys.stderr)
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'{desc} 失败 (rc={r.returncode})\n'
                           f'--- stdout ---\n{r.stdout[-3000:]}\n'
                           f'--- stderr ---\n{r.stderr[-3000:]}')


def count_mol2_topo(mol2: str) -> tuple[int, int]:
    """统计 mol2 的键/角数（用于多分子拓扑期望校验）。

    直接解析 mol2 文本（parmed 对 mol2 返回 ResidueTemplate，只含 bonds 不含 angles）。
    antechamber 生成的 mol2 通常无 @<TRIPOS>ANGLE 段（tleap 从连接性推断角），
    此时用键连接图推断：角数 = Σ 原子度 d × (d-1) / 2（与 tleap 规则一致）。
    """
    nb = na = 0
    section = ''
    bonds: list[tuple[int, int]] = []
    has_angle = False
    with open(mol2) as f:
        for ln in f:
            s = ln.strip()
            if s.startswith('@<TRIPOS>'):
                section = s
                if section == '@<TRIPOS>ANGLE':
                    has_angle = True
                continue
            if not s or s.startswith('#'):
                continue
            if section == '@<TRIPOS>BOND':
                nb += 1
                p = ln.split()
                if len(p) >= 4:
                    bonds.append((int(p[1]), int(p[2])))
            elif section == '@<TRIPOS>ANGLE':
                na += 1
    if not has_angle and bonds:
        deg: dict[int, int] = {}
        for a, b in bonds:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        na = sum(d * (d - 1) // 2 for d in deg.values())
    return nb, na


def _leaprc(forcefield: str) -> str:
    return 'leaprc.gaff2' if forcefield == 'gaff2' else 'leaprc.gaff'


# ---------------------------------------------------------------- 单分子

def build_system_single(cfg: Config, workdir: str, mol2: str, frcmod: str) -> dict:
    """单分子 tleap：mol2 + frcmod → prmtop/inpcrd。返回路径。"""
    name = cfg.name
    prmtop = os.path.join(workdir, name + '.prmtop')
    inpcrd = os.path.join(workdir, name + '.inpcrd')

    tleap_in = os.path.join(workdir, 'tleap.in')
    with open(tleap_in, 'w') as f:
        f.write(f'source {_leaprc(cfg.forcefield)}\n'
                f'mol = loadmol2 {os.path.basename(mol2)}\n'
                f'loadamberparams {os.path.basename(frcmod)}\n'
                f'saveamberparm mol {os.path.basename(prmtop)} {os.path.basename(inpcrd)}\n'
                f'quit\n')

    print(f'== 体系层: tleap (leaprc={_leaprc(cfg.forcefield)}) ==')
    _run(['tleap', '-f', os.path.basename(tleap_in)], workdir, 'tleap')
    print(f'  tleap → {os.path.relpath(prmtop)} / {os.path.relpath(inpcrd)}')
    return {'prmtop': prmtop, 'inpcrd': inpcrd}


# ---------------------------------------------------------------- 多分子

def _write_packed_pdb(cfg: Config, mol2_list: list[str], blocks: list[list[list[tuple]]],
                      out_pdb: str) -> None:
    """packed 坐标块 + 各类型 mol2 原子名 → 合并 PDB（残基名来自 mol2，TER 分隔）。

    blocks[t][k] = [(elem, x, y, z), ...]（第 t 类型第 k 个分子的坐标）
    """
    tmpl = []
    for mol2 in mol2_list:
        s = pmd.load_file(mol2)
        tmpl.append([(a.name, _element_symbol(a)) for a in s.atoms])

    lines = []
    serial = 0
    resseq = 0
    for t, block in enumerate(blocks):
        s = pmd.load_file(mol2_list[t])
        resname = _resname_of(s)
        for k, mol_coords in enumerate(block):
            resseq += 1
            for j in range(len(mol_coords)):
                serial += 1
                name, elem = tmpl[t][j]
                x, y, z = mol_coords[j][1], mol_coords[j][2], mol_coords[j][3]
                lines.append(_pdb_atom_line(serial, name, resname, resseq, x, y, z, elem))
            lines.append('TER')   # 断开分子间键合
    lines.append('END')
    with open(out_pdb, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _resname_of(s) -> str:
    """parmed 读 mol2 可能返回 ResidueTemplate（无 residues 属性），取 .name。"""
    if hasattr(s, 'residues') and s.residues:
        return s.residues[0].name
    name = getattr(s, 'name', None)
    if name:
        return str(name).strip()
    return 'MOL'


def _pdb_atom_line(serial: int, name: str, resname: str, resseq: int,
                   x: float, y: float, z: float, elem: str) -> str:
    """标准 PDB ATOM 行（tleap 可解析的列宽）。"""
    return (f'ATOM  {serial:5d} {name:>4s} {resname:>3s} A{resseq:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {elem:>2s}')


def _element_symbol(atom) -> str:
    """从原子名推元素符号（antechamber mol2 原子名 = 元素+序号，如 C1/Cl1/H1）。

    不用 parmed atomic_number：读 mol2（无质量字段）时 parmed 会把 Cl 的
    atomic_number 错识别成 6(C)，导致 PDB 元素列写成 C。
    """
    name = atom.name
    if len(name) >= 2 and name[1].islower():
        return name[:2].capitalize()
    return name[0].upper()


def _mol2_to_prepi(mol2: str, prepi: str, cwd: str) -> None:
    """mol2 → prepi（保留电荷）。teLeap 的 loadpdb 只认 prep/lib 模板，
    loadmol2 单元不注册为残基模板（AmberTools 24+ 行为差异）。"""
    _run(['antechamber', '-i', os.path.basename(mol2), '-fi', 'mol2',
          '-o', os.path.basename(prepi), '-fo', 'prepi'],
         cwd, f'antechamber({os.path.basename(prepi)})')


def build_system_multi(cfg: Config, workdir: str, mol2_list: list[str],
                       frcmod_list: list[str], blocks: list[list[list[tuple]]]) -> dict:
    """多分子 tleap：合并 PDB + loadamberprep 模板 → 体系 prmtop/inpcrd。"""
    packed_pdb = os.path.join(workdir, 'packed.pdb')
    _write_packed_pdb(cfg, mol2_list, blocks, packed_pdb)
    print(f'  合并 PDB: {os.path.relpath(packed_pdb)} '
          f'({sum(len(b) for b in blocks)} 分子)')

    # prepi 模板（teLeap loadpdb 需要）
    prepi_list = []
    for mol2 in mol2_list:
        prepi = os.path.splitext(mol2)[0] + '.prepi'
        _mol2_to_prepi(mol2, prepi, workdir)
        prepi_list.append(prepi)

    name = cfg.name
    prmtop = os.path.join(workdir, name + '.prmtop')
    inpcrd = os.path.join(workdir, name + '.inpcrd')

    tleap_in = os.path.join(workdir, 'tleap.in')
    with open(tleap_in, 'w') as f:
        f.write(f'source {_leaprc(cfg.forcefield)}\n')
        for prepi in prepi_list:
            f.write(f'loadamberprep {os.path.basename(prepi)}\n')
        for x in frcmod_list:
            f.write(f'loadamberparams {os.path.basename(x)}\n')
        f.write(f'sys = loadpdb {os.path.basename(packed_pdb)}\n')
        f.write(f'saveamberparm sys {os.path.basename(prmtop)} {os.path.basename(inpcrd)}\n')
        f.write('quit\n')

    print(f'== 体系层: tleap loadpdb (leaprc={_leaprc(cfg.forcefield)}) ==')
    _run(['tleap', '-f', os.path.basename(tleap_in)], workdir, 'tleap')
    print(f'  tleap → {os.path.relpath(prmtop)} / {os.path.relpath(inpcrd)}')
    return {'prmtop': prmtop, 'inpcrd': inpcrd, 'pdb': packed_pdb}
