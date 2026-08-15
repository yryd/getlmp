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
                      out_pdb: str, extras: list[dict] | None = None) -> None:
    """packed 坐标块 + 各类型原子名 → 合并 PDB（残基名来自 mol2，TER 分隔）。

    blocks[t][k] = [(elem, x, y, z), ...]（第 t 类型第 k 个分子的坐标）。
    extras: 溶质之外的模板类型（水/离子），与 blocks 后半段一一对应：
      {'kind': 'water', 'model': 'tip3p'}      → 残基 WAT（原子名/元素来自 off 模板）
      {'kind': 'ion', 'ion': 'Na+'}            → 残基/原子名来自 atomic_ions.lib
    """
    tmpl = []
    for mol2 in mol2_list:
        s = pmd.load_file(mol2)
        tmpl.append([(a.name, _element_symbol(a)) for a in s.atoms])

    # extras 模板（水/离子）：原子名 + 元素 + 残基名（PDB 必须与 tleap 库模板一致）
    ext_tmpl: list[dict] = []
    for ex in (extras or []):
        kind = ex['kind']
        if kind == 'water':
            from solvent_templates import water_template_xyz
            rows = water_template_xyz(ex['model'])
            ext_tmpl.append({'resname': 'WAT',
                             'atoms': [(n, e) for n, e, *_ in rows]})
        elif kind == 'ion':
            from solvent_templates import ion_type
            t = ion_type(ex['ion'])
            ext_tmpl.append({'resname': t['resname'],
                             'atoms': [(t['atomname'], t['elem'])]})
        else:
            raise RuntimeError(f'未知模板类型 kind={kind!r}（water/ion）')

    lines = []
    serial = 0
    resseq = 0
    for t, block in enumerate(blocks):
        if t < len(mol2_list):
            s = pmd.load_file(mol2_list[t])
            resname = _resname_of(s)
            atoms = tmpl[t]
        else:
            et = ext_tmpl[t - len(mol2_list)]
            resname = et['resname']
            atoms = et['atoms']
        if len(resname) > 3:
            raise RuntimeError(
                f'PDB 残基名 {resname!r} 超 3 字符（tleap loadpdb 列宽限制）；'
                f'该离子模板不支持（atomic_ions.lib 中残基名 {resname!r}）')
        for k, mol_coords in enumerate(block):
            resseq += 1
            for j in range(len(mol_coords)):
                serial += 1
                name, elem = atoms[j]
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


# 经典 TIP4P 水（AmberTools 无 leaprc.water.tip4p，getlmp 内置）：
# 残基 TP4（solvents.lib 单残基模板，经典 TIP4P 电荷 0.52/-1.04、M 距 O 0.15 Å）
# + 参数 frcmod.tip4p；离子参数沿用 JC_TIP3P
# 近似（TIP4P 无专门离子参数集，与一期 TIP4P/水 3-site 方案一致）。
# 注意：勿用 TIP4PBOX（216 水盒子模板）——loadpdb 会整盒复制导致原子数暴增。
# addAtomTypes 离子列表照抄 leaprc.water.tip4pew（loadpdb 需类型→元素映射）。
#
# 注意：AmberTools 的 frcmod.tip4p 只有 OW-EP 键，没有 OW-HW / HW-HW 键参数
# （TIP4P-Ew/OPC 等的 frcmod 均有 OW-HW 553.0/0.9572 与 HW-HW 553.0/1.5136）
# → 需 getlmp 内置补充 frcmod，否则 tleap 报 "Could not find bond parameter
# OW-HW / HW-HW" 且 saveamberparm 失败（prmtop 不生成）。
# HW-HW 是 TIP4PBOX 模板 connectivity 自带键（H1-H2），ParmEd 修复阶段会删。
_TIP4P_OWHW_FRCMOD = '''\
This is getlmp additional parameters for classic TIP4P water (OW-HW/HW-HW bonds; AmberTools frcmod.tip4p missing them)
MASS

BOND
OW-HW   553.000   0.9572
HW-HW   553.000   1.5136

ANGLE

DIHE

NONBON

'''

_TIP4P_LEAPRC = '''\
# getlmp 内置：经典 TIP4P（AmberTools 无 leaprc.water.tip4p）
addAtomTypes {
\t{ "OW"   "O"  "sp3" }
\t{ "HW"   "H"  "sp3" }
\t{ "EP"   ""   "sp3" }
\t{ "F-"   "F"  "sp3" }
\t{ "Cl-"  "Cl" "sp3" }
\t{ "Br-"  "Br" "sp3" }
\t{ "I-"   "I"  "sp3" }
\t{ "Li+"  "Li" "sp3" }
\t{ "Na+"  "Na" "sp3" }
\t{ "K+"   "K"  "sp3" }
\t{ "Rb+"  "Rb" "sp3" }
\t{ "Cs+"  "Cs" "sp3" }
\t{ "Mg+"  "Mg" "sp3" }
\t{ "Tl+"  "Tl" "sp3" }
\t{ "Cu+"  "Cu" "sp3" }
\t{ "Ag+"  "Ag" "sp3" }
\t{ "Be2+" "Be" "sp3" }
\t{ "Cu2+" "Cu" "sp3" }
\t{ "Ni2+" "Ni" "sp3" }
\t{ "Pt2+" "Pt" "sp3" }
\t{ "Zn2+" "Zn" "sp3" }
\t{ "Co2+" "Co" "sp3" }
\t{ "Pd2+" "Pd" "sp3" }
\t{ "Ag2+" "Ag" "sp3" }
\t{ "Cr2+" "Cr" "sp3" }
\t{ "Fe2+" "Fe" "sp3" }
\t{ "Mg2+" "Mg" "sp3" }
\t{ "V2+"  "V"  "sp3" }
\t{ "Mn2+" "Mn" "sp3" }
\t{ "Hg2+" "Hg" "sp3" }
\t{ "Cd2+" "Cd" "sp3" }
\t{ "Yb2+" "Yb" "sp3" }
\t{ "Ca2+" "Ca" "sp3" }
\t{ "Sn2+" "Sn" "sp3" }
\t{ "Pb2+" "Pb" "sp3" }
\t{ "Eu2+" "Eu" "sp3" }
\t{ "Sr2+" "Sr" "sp3" }
\t{ "Sm2+" "Sm" "sp3" }
\t{ "Ba2+" "Ba" "sp3" }
\t{ "Ra2+" "Ra" "sp3" }
\t{ "Al3+" "Al" "sp3" }
\t{ "Fe3+" "Fe" "sp3" }
\t{ "Cr3+" "Cr" "sp3" }
\t{ "In3+" "In" "sp3" }
\t{ "Tl3+" "Tl" "sp3" }
\t{ "Y3+"  "Y"  "sp3" }
\t{ "La3+" "La" "sp3" }
\t{ "Ce3+" "Ce" "sp3" }
\t{ "Pr3+" "Pr" "sp3" }
\t{ "Nd3+" "Nd" "sp3" }
\t{ "Sm3+" "Sm" "sp3" }
\t{ "Eu3+" "Eu" "sp3" }
\t{ "Gd3+" "Gd" "sp3" }
\t{ "Tb3+" "Tb" "sp3" }
\t{ "Dy3+" "Dy" "sp3" }
\t{ "Er3+" "Er" "sp3" }
\t{ "Tm3+" "Tm" "sp3" }
\t{ "Lu3+" "Lu" "sp3" }
\t{ "Hf4+" "Hf" "sp3" }
\t{ "Zr4+" "Zr" "sp3" }
\t{ "Ce4+" "Ce" "sp3" }
\t{ "U4+"  "U"  "sp3" }
\t{ "Pu4+" "Pu" "sp3" }
\t{ "Th4+" "Th" "sp3" }
}
loadOff atomic_ions.lib
loadOff solvents.lib
WAT = TP4
loadAmberParams frcmod.tip4p
loadAmberParams frcmod.tip4p_owhw
loadAmberParams frcmod.ionsjc_tip3p
loadAmberParams frcmod.ions234lm_1264_tip3p
'''


def _write_tip4p_leaprc(path: str) -> None:
    """写出 getlmp 内置经典 TIP4P leaprc + 补充 frcmod（AmberTools 无现成）。"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(_TIP4P_LEAPRC)
    frcmod = os.path.join(os.path.dirname(path), 'frcmod.tip4p_owhw')
    with open(frcmod, 'w', encoding='utf-8') as f:
        f.write(_TIP4P_OWHW_FRCMOD)


def _fix_water_topology(prmtop: str, inpcrd: str, model: str) -> dict:
    """修复水分子拓扑（AmberTools 水模板无 connect/angle，loadpdb 需校正）。

    背景：solvents.lib 的 TP3/SPC/OPC3 水模板没有 connect0（分子内键）和
    angle 定义——Amber 水靠 SETTLE 约束，sander/pmemd 不需要显式键角参数。
    tleap loadpdb 对无模板键的残基走自动键合：
      - O-H1 / O-H2 键 ✓（正确）
      - H1-H2 误建键 ✗（H-H 距离 1.51 Å 落入自动键合阈值）
      - 角 0 ✗（tleap 角只来自模板定义，不从键图推导）
    这里用 ParmEd 校正：删 H1-H2 键 + 确保 O-H1/O-H2 键 + 补 H1-O-H2 角
    （角类型 HW-OW-HW，参数取水模型表的标准值；LAMMPS 需要显式角参数，
    刚性水模拟可再用 fix rigid/settle，角参数被忽略）。
    返回统计 {'bonds_removed', 'bonds_added', 'angles_added'}。
    """
    from parmed.topologyobjects import Angle, AngleType, Bond
    from solvent_templates import water_model

    wm = water_model(model)
    s = pmd.load_file(prmtop, inpcrd)
    stat = {'bonds_removed': 0, 'bonds_added': 0, 'angles_added': 0}

    # 1) 删 WAT 残基内误建键：
    #    - H1-H2（自动键合误建，H-H 距离 1.51 Å 落入阈值）
    #    - O-EPW（4-site 虚拟位点键；LAMMPS 隐式 M 方案不需要 EP 原子/键，
    #      导出层会丢弃 EPW 并把电荷合并到 O）
    for b in list(s.bonds):
        a1, a2 = b.atom1, b.atom2
        if (a1.residue.name == 'WAT' and a2.residue.name == 'WAT'
                and {a1.name, a2.name} in ({'H1', 'H2'}, {'O', 'EPW'})):
            s.bonds.remove(b)
            stat['bonds_removed'] += 1

    # 2) 确保 O-H1/O-H2 键存在（自动键合通常已建；缺失时补，用现有 OW-HW 类型）
    ow_type = None
    for b in s.bonds:
        if b.atom1.residue.name == 'WAT' and {b.atom1.name, b.atom2.name} == {'O', 'H1'}:
            ow_type = b.type
            break
    for res in s.residues:
        if res.name != 'WAT':
            continue
        atoms = {a.name: a for a in res.atoms}
        o = atoms['O']
        for hn in ('H1', 'H2'):
            h = atoms[hn]
            if not any((b.atom1 is o and b.atom2 is h)
                       or (b.atom1 is h and b.atom2 is o) for b in s.bonds):
                s.bonds.append(Bond(o, h, ow_type))
                stat['bonds_added'] += 1

    # 3) 补 H1-O-H2 角（HW-OW-HW；k/θeq 取水模型表标准值）
    atype = AngleType(wm['angle_k'], wm['angle_theta'])
    s.angle_types.append(atype)
    for res in s.residues:
        if res.name != 'WAT':
            continue
        atoms = {a.name: a for a in res.atoms}
        o, h1, h2 = atoms['O'], atoms['H1'], atoms['H2']
        if not any((a.atom1 is h1 and a.atom2 is o and a.atom3 is h2)
                   or (a.atom1 is h2 and a.atom2 is o and a.atom3 is h1)
                   for a in s.angles):
            s.angles.append(Angle(h1, o, h2, atype))
            stat['angles_added'] += 1

    # parmed save 不覆盖已存在文件 → 写临时再原子替换
    tmp_p, tmp_i = prmtop + '.fix', inpcrd + '.fix'
    s.save(tmp_p, format='amber')
    s.save(tmp_i, format='rst7')
    os.replace(tmp_p, prmtop)
    os.replace(tmp_i, inpcrd)
    return stat


def build_system_multi(cfg: Config, workdir: str, mol2_list: list[str],
                       frcmod_list: list[str], blocks: list[list[list[tuple]]],
                       extras: list[dict] | None = None) -> dict:
    """多分子 tleap：合并 PDB + loadamberprep 模板 → 体系 prmtop/inpcrd。

    extras: 水/离子模板类型（见 _write_packed_pdb）；水存在时 tleap 额外
    `source leaprc.water.{model}`（加载水 + 离子 LJ 参数，与 yaml 选择的
    水模型配套）。无 extras 时行为与旧版一致。
    """
    packed_pdb = os.path.join(workdir, 'packed.pdb')
    _write_packed_pdb(cfg, mol2_list, blocks, packed_pdb, extras=extras)
    print(f'  合并 PDB: {os.path.relpath(packed_pdb)} '
          f'({sum(len(b) for b in blocks)} 分子'
          + (f' + {len(extras)} 模板类型' if extras else '') + ')')

    # prepi 模板（teLeap loadpdb 需要）
    prepi_list = []
    for mol2 in mol2_list:
        prepi = os.path.splitext(mol2)[0] + '.prepi'
        _mol2_to_prepi(mol2, prepi, workdir)
        prepi_list.append(prepi)

    # 水模型 leaprc（若含水的模板类型；水模型决定配水离子参数 frcmod.ions*；
    # 经典 TIP4P 无现成 leaprc → getlmp 内置生成）
    water_leaprc = None
    water_model_name = None
    for ex in (extras or []):
        if ex['kind'] == 'water':
            from solvent_templates import water_model
            water_model_name = ex['model']
            water_leaprc = water_model(water_model_name)['leaprc']
            break

    name = cfg.name
    prmtop = os.path.join(workdir, name + '.prmtop')
    inpcrd = os.path.join(workdir, name + '.inpcrd')

    tleap_in = os.path.join(workdir, 'tleap.in')
    with open(tleap_in, 'w') as f:
        f.write(f'source {_leaprc(cfg.forcefield)}\n')
        if water_leaprc:
            wm = water_model(water_model_name)
            if wm.get('custom_leaprc'):
                custom = os.path.join(workdir, 'leaprc.water.tip4p')
                _write_tip4p_leaprc(custom)
                f.write(f'source {os.path.basename(custom)}\n')
            else:
                f.write(f'source {water_leaprc}\n')
        for prepi in prepi_list:
            f.write(f'loadamberprep {os.path.basename(prepi)}\n')
        for x in frcmod_list:
            f.write(f'loadamberparams {os.path.basename(x)}\n')
        f.write(f'sys = loadpdb {os.path.basename(packed_pdb)}\n')
        f.write(f'saveamberparm sys {os.path.basename(prmtop)} {os.path.basename(inpcrd)}\n')
        f.write('quit\n')

    print(f'== 体系层: tleap loadpdb (leaprc={_leaprc(cfg.forcefield)}'
          + (f' + {water_leaprc}'
             + (' [getlmp 内置]' if water_leaprc and water_model(water_model_name).get('custom_leaprc') else '')
             if water_leaprc else '') + ') ==')
    _run(['tleap', '-f', os.path.basename(tleap_in)], workdir, 'tleap')
    print(f'  tleap → {os.path.relpath(prmtop)} / {os.path.relpath(inpcrd)}')

    # 水拓扑修复（AmberTools 水模板无 connect/angle，见 _fix_water_topology）
    water_fix = None
    if water_leaprc:
        model = next(ex['model'] for ex in (extras or []) if ex['kind'] == 'water')
        water_fix = _fix_water_topology(prmtop, inpcrd, model)
        wm = water_model(model)
        print(f'  水拓扑修复: 删误建键 {water_fix["bonds_removed"]}, '
              f'补 O-H 键 {water_fix["bonds_added"]}, '
              f'补 H-O-H 角 {water_fix["angles_added"]}'
              f'（角参数 HW-OW-HW k={wm["angle_k"]:g} '
              f'θeq={wm["angle_theta"]:g}°）')
    return {'prmtop': prmtop, 'inpcrd': inpcrd, 'pdb': packed_pdb,
            'water_fix': water_fix}
