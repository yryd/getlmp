#!/usr/bin/env python
"""体系层·packmol：多分子装盒。

输入：各分子类型的 mol2（含坐标/原子名）+ packmol 配置
输出：packed.xyz（按 structure 顺序排列的全部坐标）
"""
from __future__ import annotations

import os
import subprocess
import sys

import parmed as pmd

from config import Config, PackmolCfg


def mol2_to_xyz(mol2: str, out_xyz: str) -> int:
    """mol2 → packmol 输入 xyz（元素 + 坐标，原子顺序与 mol2 一致）。返回原子数。"""
    s = pmd.load_file(mol2)
    lines = [str(len(s.atoms)), 'from getlmp (mol2)']
    for a in s.atoms:
        elem = _element_symbol(a)
        lines.append(f'{elem:2s} {a.xx:12.5f} {a.xy:12.5f} {a.xz:12.5f}')
    with open(out_xyz, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return len(s.atoms)


def sdf_to_xyz(sdf: str, out_xyz: str) -> int:
    """SDF → packmol 输入 xyz（元素 + 坐标）。用于 ReaxFF 分支（无 mol2）。返回原子数。"""
    from rdkit import Chem
    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        raise RuntimeError(f'RDKit 无法读取 SDF: {sdf}')
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), 'from getlmp (sdf)']
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f'{a.GetSymbol():2s} {p.x:12.5f} {p.y:12.5f} {p.z:12.5f}')
    with open(out_xyz, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return mol.GetNumAtoms()


def _element_symbol(atom) -> str:
    z = getattr(atom, 'atomic_number', 0)
    if z and z > 0:
        return _ELEM[z] if z < len(_ELEM) else atom.name[:1].upper()
    # 兜底：从原子名推元素（antechamber mol2 原子名通常为 元素+序号）
    name = atom.name
    if len(name) >= 2 and name[1].islower():
        return name[:2].capitalize()
    return name[0].upper()


_ELEM = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']


def write_inp(cfg: Config, struct_xyz: list[str], n_atoms: list[int],
              inp_path: str) -> None:
    """生成 packmol inp（阶段 2：bulk 预设，每类型一个 structure 块 inside box）。"""
    pm: PackmolCfg = cfg.packmol
    xlo, ylo, zlo, xhi, yhi, zhi = pm.box
    lines = [
        f'tolerance {pm.tolerance:.2f}',
        f'seed {pm.seed}',
        'filetype xyz',
        'output packed.xyz',
        '',
    ]
    for i, mc in enumerate(cfg.molecules):
        lines += [
            f'structure {os.path.basename(struct_xyz[i])}',
            f'  number {mc.count}',
            f'  inside box {xlo:.3f} {ylo:.3f} {zlo:.3f} {xhi:.3f} {yhi:.3f} {zhi:.3f}',
            'end structure',
            '',
        ]
    with open(inp_path, 'w') as f:
        f.write('\n'.join(lines))


def parse_inp_structures(inp_path: str) -> list[tuple[str, int]]:
    """解析 packmol inp → [(structure 文件名, number), ...]（按出现顺序）。

    用于自定义 inp（packmol.inp_file）：structure 顺序须与 molecules 一致
    （用户保证），number 以 inp 内为准（可与 yaml count 不同，分块/守恒校验
    都用它）。
    """
    with open(inp_path, encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    structs: list[tuple[str, int]] = []
    cur: str | None = None
    for ln in lines:
        s = ln.strip()
        if s.lower().startswith('structure '):
            cur = s.split(None, 1)[1].strip()
        elif s.lower().startswith('number ') and cur is not None:
            try:
                n = int(s.split(None, 1)[1].strip())
            except ValueError:
                raise RuntimeError(f'packmol inp number 解析失败: {ln!r}')
            structs.append((cur, n))
            cur = None
    if not structs:
        raise RuntimeError(f'packmol inp 未找到 structure/number 块: {inp_path}')
    return structs


def run_packmol(inp_path: str, workdir: str) -> str:
    """跑 packmol（21.x 用法 `packmol < inp`）。

    注意：packmol 是 Fortran 程序，stdin 必须是可 seek 的真实文件
    （不能用 subprocess 的 pipe），否则报 'Illegal seek'。
    """
    print(f'  [run] packmol < {os.path.basename(inp_path)}', file=sys.stderr)
    with open(inp_path) as f_in:
        r = subprocess.run(['packmol'], stdin=f_in, cwd=workdir,
                           capture_output=True, text=True, timeout=600)
    out = (r.stdout or '') + (r.stderr or '')
    if r.returncode != 0:
        raise RuntimeError(f'packmol 失败 (rc={r.returncode})\n{out[-3000:]}')
    # packmol 成功输出 "Success!"（21.x 用感叹号；兼容旧版 'Successfully'）
    if 'Success' not in out:
        raise RuntimeError(f'packmol 未成功结束（未见 Success）\n{out[-3000:]}')
    packed = os.path.join(workdir, 'packed.xyz')
    if not os.path.exists(packed):
        raise RuntimeError(f'packmol 未生成 {packed}\n{out[-3000:]}')
    return packed


def parse_packed_xyz(packed_xyz: str, counts: list[int],
                     natom_by_name: dict[str, int],
                     type_names: list[str]) -> list[list[list[tuple]]]:
    """解析 packed.xyz → 每个分子类型一个坐标块。

    返回 list[type_index] -> list[分子拷贝] -> list[(elem, x, y, z)]。
    packmol 输出按 structure 顺序连续排列，块大小 = number × natom(type)。
    counts/type_names 一一对应（同一类型可拆多个 structure 块，如"固定 1 个 +
    自由 N 个"）；blocks 按类型首次出现顺序排列（与 molecules 顺序一致的前提：
    structure 首块顺序 = molecules 顺序）。
    """
    with open(packed_xyz) as f:
        lines = f.read().splitlines()
    total = int(lines[0].split()[0])
    expect = sum(c * natom_by_name[t] for c, t in zip(counts, type_names))
    if total != expect:
        raise RuntimeError(f'packed.xyz 原子数 {total} != 期望 {expect}')
    coords = []
    for ln in lines[2:]:
        p = ln.split()
        if len(p) < 4:
            continue
        coords.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    if len(coords) != total:
        raise RuntimeError(f'packed.xyz 实际坐标行 {len(coords)} != 头部 {total}')

    order: list[str] = []
    for t in type_names:
        if t not in order:
            order.append(t)
    blocks = {t: [] for t in order}
    idx = 0
    for c, t in zip(counts, type_names):
        n = natom_by_name[t]
        for _ in range(c):
            blocks[t].append(coords[idx:idx + n])
            idx += n
    return [blocks[t] for t in order]
