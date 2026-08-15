#!/usr/bin/env python
"""导出层：prmtop/inpcrd → data.lmp（自写转换器）+ 格式/守恒校验。"""
from __future__ import annotations

import os

from check_lammps_data import parse_data
from config import Config
from prmtop_to_lammps import prmtop_to_lammps


def build_data(cfg: Config, workdir: str, prmtop: str, inpcrd: str,
               box: list | None = None) -> dict:
    """prmtop/inpcrd → data.lmp。返回 (导出 info, 自检 result)。

    box: 显式盒 [xlo,xhi,ylo,yhi,zlo,zhi]（多分子由 packmol 决定；单分子用 buffer）
    """
    out_lmp = os.path.join(workdir, cfg.output)
    print('== 导出层: prmtop → data.lmp ==')
    info = prmtop_to_lammps(prmtop, inpcrd, out_lmp, buffer=cfg.buffer, box=box)
    check = parse_data(out_lmp)
    print(f'  data.lmp: {info["natom"]} atoms / {info["nbond"]} bonds / '
          f'{info["nangle"]} angles / {info["ndihedral"]} dihedrals / '
          f'{info["nimproper"]} impropers / charge={info["total_charge"]:.6f}')
    return {'info': info, 'check': check, 'data_lmp': out_lmp}


def validate_export(cfg: Config, result: dict, natom_input: int,
                    expected_topo: dict | None = None) -> tuple[bool, list[str]]:
    """守恒与一致性校验：电荷、原子数、段计数、格式自检。

    expected_topo: {'bonds': N, 'angles': M} 期望拓扑数（多分子：Σ count×单分子键/角数），
                   用于验证分子键合完整（如水分子 O-H 键齐全）。
    """
    info = result['info']
    check = result['check']
    msgs: list[str] = []
    ok = True

    # 1) 电荷守恒（容差 0.01）
    dq = abs(info['total_charge'] - cfg.net_charge)
    if dq <= 0.01:
        msgs.append(f'电荷守恒: sum={info["total_charge"]:.6f} 净电荷={cfg.net_charge} 通过')
    else:
        ok = False
        msgs.append(f'电荷守恒: sum={info["total_charge"]:.6f} != 净电荷={cfg.net_charge} 失败')

    # 2) 原子数守恒（分子层含氢原子数；4-site 水 EP 虚拟位点在导出层已
    #    丢弃并合并电荷 → 期望 = 分子层原子数 - EP 数）
    ep_removed = info.get('ep_removed', 0)
    natom_expect = natom_input - ep_removed
    ep_note = f'（-{ep_removed} EP 虚拟位点）' if ep_removed else ''
    if info['natom'] == natom_expect:
        msgs.append(f'原子数守恒: {info["natom"]} == 分子层 {natom_input}'
                    f'{ep_note} 通过')
    else:
        ok = False
        msgs.append(f'原子数不一致: data={info["natom"]} 分子层={natom_input}'
                    f'{ep_note} 失败')

    # 3) 段计数（头部 vs 实际）
    for key in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers'):
        h = check['counts'].get(key, 0)
        a = check.get(key + '_actual', 0)
        if h != a:
            ok = False
            msgs.append(f'段计数不一致 {key}: 头部={h} 实际={a} 失败')
    msgs.append('段计数: 头部与文件记录一致 通过')

    # 4) 拓扑期望（多分子键合完整：键/角数 == Σ count×单分子拓扑）
    if expected_topo:
        topo_key = {'bonds': 'nbond', 'angles': 'nangle'}
        for key, exp in expected_topo.items():
            act = info[topo_key[key]]
            if act == exp:
                msgs.append(f'拓扑期望 {key}: {act} == {exp} 通过')
            else:
                ok = False
                msgs.append(f'拓扑期望 {key}: data={act} != 期望 {exp} 失败（分子键合可能不完整）')

    # 5) 类型计数
    for key in ('atom_types', 'bond_types', 'angle_types', 'dihedral_types', 'improper_types'):
        h = check['counts'].get(key, 0)
        a = check.get(key + '_actual', 0)
        if h != a:
            ok = False
            msgs.append(f'类型计数不一致 {key}: 头部={h} 实际={a} 失败')
    msgs.append('类型计数: 头部与 Coeffs 段一致 通过')

    return ok, msgs
