#!/usr/bin/env python
"""体系层·packmol：多分子装盒。

输入：各分子类型的 mol2（含坐标/原子名）+ packmol 配置
输出：packed.xyz（按 structure 顺序排列的全部坐标）
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
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
    """从原子名推元素符号（antechamber mol2 原子名 = 元素+序号，如 C1/Cl1/H1/C9）。

    不用 parmed atomic_number：读 mol2（无质量字段）时 parmed 会把 Cl 的
    atomic_number 错识别成 6(C)，导致元素列/密度盒子算错。
    """
    name = atom.name
    if len(name) >= 2 and name[1].islower():
        return name[:2].capitalize()
    return name[0].upper()


def write_inp(cfg: Config, structs: list[tuple[str, str, int]], inp_path: str) -> None:
    """生成 packmol inp（阶段 2：bulk 预设，每类型一个 structure 块 inside box）。

    structs: [(类型名, xyz 文件路径, 拷贝数), ...]，顺序 = packmol 输出顺序
    （溶质 + 水 + 离子；与合并 PDB / 拓扑期望校验一致）。

    inside box 缩进 padding = tolerance/2：packmol 只做欧氏距离检查（不感知
    周期边界），分子可贴盒边界导致 PBC 下相邻镜像重叠（曾致 data.lmp 能量
    爆炸 ~1e8 kcal/mol）。缩进后所有原子距边界 ≥ tol/2，任意跨边界分子对的
    周期最小镜像距离 ≥ tol/2 + tol/2 = tolerance（数学保证无 PBC 重叠）。
    """
    pm: PackmolCfg = cfg.packmol
    xlo, ylo, zlo, xhi, yhi, zhi = pm.box
    pad = pm.tolerance / 2.0
    lines = [
        f'tolerance {pm.tolerance:.2f}',
        f'nloop0 {pm.nloop0:d}',
        f'seed {pm.seed}',
        'filetype xyz',
        'output packed.xyz',
        '',
    ]
    for name, xyz_path, count in structs:
        lines += [
            f'structure {os.path.basename(xyz_path)}',
            f'  number {count}',
            f'  inside box {xlo + pad:.3f} {ylo + pad:.3f} {zlo + pad:.3f} '
            f'{xhi - pad:.3f} {yhi - pad:.3f} {zhi - pad:.3f}',
            'end structure',
            '',
        ]
    with open(inp_path, 'w') as f:
        f.write('\n'.join(lines))


AVOGADRO = 6.02214076e23
# 相对原子质量（标准原子量，与 tleap mass 近似一致，仅用于密度盒子估算）
ELEMENT_MASS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'Na': 22.990, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'K': 39.098,
    'Ca': 40.078, 'Br': 79.904, 'I': 126.904,
}


def mol2_mass(mol2: str) -> float:
    """从 mol2 读单分子摩尔质量（g/mol）——按元素符号查 ELEMENT_MASS。

    注意不能用 parmed atom.mass：mol2 无质量字段，parmed 按 ResidueTemplate
    读出的 mass 为 0（曾导致密度盒子算出 0）。
    """
    s = pmd.load_file(mol2)
    total = 0.0
    for a in s.atoms:
        elem = _element_symbol(a)
        m = ELEMENT_MASS.get(elem)
        if m is None:
            raise RuntimeError(f'密度盒子无法计算：mol2 中未知元素 {elem!r}'
                               f'（ELEMENT_MASS 未收录，见 src/packmol_layer.py）')
        total += m
    return total


def elements_mass(elements: list[str]) -> float:
    """按元素组成估算摩尔质量（g/mol）——ReaxFF 分支用（无 mol2 质量表）。"""
    total = 0.0
    for e in elements:
        m = ELEMENT_MASS.get(e)
        if m is None:
            raise RuntimeError(f'密度盒子无法计算：未知元素 {e!r}（ELEMENT_MASS 未收录）')
        total += m
    return total


def density_box(total_mass_g: float, density: float) -> list:
    """按 总质量(g)/密度(g/cm³) → 立方盒 [0,0,0,L,L,L]（Å）。

    1 g/cm³ = 1e-24 g/Å³（1 cm³ = 1e24 Å³）。
    与 packmol 的 density 关键字等价（packmol 默认立方盒）。
    """
    if density <= 0:
        raise ValueError(f'density 需为正数，当前 {density}')
    vol_a3 = total_mass_g / (density * 1e-24)
    L = vol_a3 ** (1.0 / 3.0)
    return [0.0, 0.0, 0.0, L, L, L]


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


# ---------------------------------------------------------------- PBC 重叠修复
# packmol 只做欧氏距离检查：分子可紧贴盒边界，PBC 下与相邻镜像重叠；
# 且球排除近似对大分子不精确，溶质-溶剂可贴脸（< tolerance）。
# 这里做周期最小镜像距离检查，把重叠的溶剂分子随机平移重排到无重叠位置。


def refine_overlap(blocks: list, type_names: list[str], box: list,
                   tol: float = 2.0, seed: int = 2026,
                   moveable_names: set | None = None) -> list:
    """PBC-aware 重叠后处理：消除 packmol 装盒残留的重叠/贴脸接触。

    blocks: list[type_index] -> list[分子拷贝] -> list[(elem, x, y, z)]
    type_names: 与 blocks 对应的类型名（'water' / 'ion_*' / 溶质 name）
    box: packmol box 语义 [xlo,ylo,zlo,xhi,yhi,zhi]
    tol: 允许的最小非键原子距离（Å），默认 2.0
    moveable_names: 可重排的类型名集合（默认 None → 全部可重排）。
      建议传溶剂类型（water/ion_*），溶质固定。

    检查所有分子对的周期最小镜像原子距离，< tol 的分子中可动者随机平移
    （保持取向）重排，直到全局无重叠或达到轮次上限。返回修正后的 blocks
    （顺序与原子数不变，仅坐标可能改变）。
    """
    if moveable_names is None:
        moveable_names = set(type_names)

    # 展平为分子列表：(type_name, 原子坐标 (n,3), 半径 r_max)
    mols: list[tuple[str, np.ndarray, float]] = []
    for ti, blk in enumerate(blocks):
        tname = type_names[ti]
        for mol in blk:
            arr = np.array([[a[1], a[2], a[3]] for a in mol], dtype=float)
            center = arr.mean(axis=0)
            r = float(np.max(np.linalg.norm(arr - center, axis=1))) if len(arr) else 0.0
            mols.append((tname, arr, r))
    L = np.array([box[3] - box[0], box[4] - box[1], box[5] - box[2]], dtype=float)
    rng = np.random.default_rng(seed)

    for _round in range(60):
        pairs = _overlap_pairs(mols, L, tol)
        if not pairs:
            break
        to_move = set()
        for i, j in pairs:
            mi, mj = mols[i][0], mols[j][0]
            fi, fj = mi not in moveable_names, mj not in moveable_names
            if fi and fj:
                continue          # 两个都固定（理论上不该发生：固定对重叠）
            if fi:
                to_move.add(j)
            elif fj:
                to_move.add(i)
            else:
                to_move.add(j)    # 都动时移后一个，保留前一个参考
        for mi in sorted(to_move):
            _relocate_molecule(mols, mi, L, tol, rng)
        if _round == 59:
            n_left = len(_overlap_pairs(mols, L, tol))
            raise RuntimeError(
                f'PBC 重叠重排未收敛（仍剩 {n_left} 对 < {tol:.1f} Å）：'
                f'体系过密，请增大盒子或减少分子数')

    # 写回 blocks（坐标顺序不变）
    idx = 0
    for ti, blk in enumerate(blocks):
        tname = type_names[ti]
        for k, mol in enumerate(blk):
            arr = mols[idx][1]
            for a_i, (elem, _x, _y, _z) in enumerate(mol):
                blk[k][a_i] = (elem, float(arr[a_i, 0]), float(arr[a_i, 1]),
                               float(arr[a_i, 2]))
            idx += 1
    return blocks


def _pbc_delta(d: np.ndarray, L: np.ndarray) -> np.ndarray:
    """周期最小镜像差值：Δ - L*round(Δ/L)。"""
    return d - L * np.round(d / L)


def _overlap_pairs(mols: list, L: np.ndarray,
                   tol: float) -> list[tuple[int, int]]:
    """找所有最小镜像原子距离 < tol 的分子对（跳过同分子）。"""
    M = len(mols)
    if M < 2:
        return []
    centers = np.array([m[1].mean(axis=0) for m in mols])       # (M,3)
    radii = np.array([m[2] for m in mols])                       # (M,)
    pairs: list[tuple[int, int]] = []
    # 质心预筛（分块避免超大矩阵）：候选 = 质心最小镜像距离 < tol + 2*max(r_i,r_j)
    max_r = radii.max() if len(radii) else 0.0
    cutoff = tol + 2.0 * max_r
    for i in range(M):
        d = centers[i + 1:] - centers[i]                         # (M-i-1,3)
        d = _pbc_delta(d, L)
        dist = np.sqrt((d * d).sum(axis=1))
        cand = np.nonzero(dist < cutoff)[0] + (i + 1)
        for j in cand:
            if _mol_pair_min_dist(mols[i][1], mols[j][1], L) < tol:
                pairs.append((i, int(j)))
    return pairs


def _mol_pair_min_dist(a: np.ndarray, b: np.ndarray,
                       L: np.ndarray) -> float:
    """两个分子间最小周期镜像原子距离。"""
    d = a[:, None, :] - b[None, :, :]                            # (na,nb,3)
    d = _pbc_delta(d, L)
    return float(np.sqrt((d * d).sum(axis=2)).min())


def _relocate_molecule(mols: list, mi: int, L: np.ndarray, tol: float,
                       rng: np.random.Generator) -> None:
    """把分子 mi 随机平移到无重叠位置（保持取向，尝试多次）。"""
    tname, arr, r = mols[mi]
    rel = arr - arr.mean(axis=0)                                 # 相对质心
    for _try in range(300):
        # 质心采样：保证分子完整在盒内（原子坐标 ∈ [0, L]）
        lo = r + 1e-6
        hi = L - r - 1e-6
        if np.any(hi <= lo):
            # 分子比盒还大（异常），退回随机平移后 wrap
            center = rng.uniform(0, L, size=3)
        else:
            center = rng.uniform(lo, hi)
        new = rel + center
        ok = True
        for j, (_tj, _aj, _rj) in enumerate(mols):
            if j == mi:
                continue
            if _mol_pair_min_dist(new, _aj, L) < tol:
                ok = False
                break
        if ok:
            mols[mi] = (tname, new, r)
            return
    # 300 次失败：接受最后一次尝试的位置（不抛出，外层轮次会再检查）
    mols[mi] = (tname, new, r)
