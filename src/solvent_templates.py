#!/usr/bin/env python
"""溶剂/离子模板：水模型与单原子离子的内置模板 + 坐标提取。

水模型：AmberTools 原生（leaprc.water.* + box off），一期 tip3p/spce/opc3，
二期预留 opc/tip4pew/tip4pd 等（supported=False 仅校验拦截）。
单分子水坐标直接从 AmberTools off 库文本解析（与 tleap 加载的库严格一致，
不依赖 parmed —— parmed 读 off 有 residue table 警告噪音）。

离子：atomic_ions.lib 动态扫描（运行时读取，参数来源与 tleap 完全一致），
单原子离子支持；多原子离子（H3O+/NH4+ 等）按需求排除。
"""
from __future__ import annotations

import os
import re
import sys

# AmberTools 数据目录（conda 环境；AMBERHOME 未设时用 sys.prefix）
def _amber_dat() -> str:
    env = os.environ.get('AMBERHOME', '')
    if env and os.path.isdir(os.path.join(env, 'dat')):
        return os.path.join(env, 'dat')
    return os.path.join(sys.prefix, 'dat')


AMBER_DAT = _amber_dat()
LIB_DIR = os.path.join(AMBER_DAT, 'leap', 'lib')
CMD_DIR = os.path.join(AMBER_DAT, 'leap', 'cmd')

# ---------------------------------------------------------------- 水模型表
# nsite: PDB 中写的原子数（3-site = 3；4-site = 4，含虚拟位点 EPW）。
# 4-site（TIP4P*/OPC）：AmberTools 的负电荷全在 EPW（O 电荷=0），
# LAMMPS 用隐式 M 方案（pair_style lj/cut/tip4p/long 运行时从 O-H 几何推导
# M 位点，data 不需要 EP 原子）→ 导出时 EP 电荷合并到 O、丢弃 EP 原子。
# om_dist: O→M 距离（LAMMPS pair_style qdist 参数；文献值，报告对照用，
# 实现时以 off 坐标实测为准）。
# custom_leaprc: AmberTools 无现成 leaprc（经典 TIP4P），getlmp 内置。
WATER_MODELS: dict[str, dict] = {
    # angle_k/angle_theta: HW-OW-HW 角参数（AmberTools 水模板无角定义，
    # SETTLE 约束设计；LAMMPS data 需要显式角参数，取各模型标准值：
    #   TIP3P: Jorgensen 100/104.52；SPC/E: Berendsen 100/109.47
    #   （frcmod.spce 键长 1.0 配套文献几何）；OPC3: 100/109.43（Izadi 2016）
    #   TIP4P*/TIP4P-D: 104.52；OPC: 103.60（LAMMPS Howto_tip4p 表）
    # 刚性水模拟可再用 fix rigid/settle，角参数被忽略）
    'tip3p': dict(leaprc='leaprc.water.tip3p', box_off='tip3pbox.off',
                  resname='WAT', nsite=3, supported=True,
                  angle_k=100.0, angle_theta=104.52,
                  desc='TIP3P（经典 3-site，GAFF 默认配水）'),
    'spce': dict(leaprc='leaprc.water.spce', box_off='spcebox.off',
                 resname='WAT', nsite=3, supported=True,
                 angle_k=100.0, angle_theta=109.47,
                 desc='SPC/E（3-site，刚性，常用）'),
    'opc3': dict(leaprc='leaprc.water.opc3', box_off='opc3box.off',
                 resname='WAT', nsite=3, supported=True,
                 angle_k=100.0, angle_theta=109.43,
                 desc='OPC3（3-site，OPC 系列简化版，模拟精度好）'),
    # ---- 二期：4-site 水（LAMMPS 隐式 M 方案，导出层 EP 合并/丢弃） ----
    'tip4p': dict(leaprc='leaprc.water.tip4p', box_off='tip4pbox.off',
                  resname='WAT', nsite=4, supported=True, has_ep=True,
                  custom_leaprc=True,
                  angle_k=100.0, angle_theta=104.52, om_dist=0.1500,
                  desc='TIP4P（经典 4-site，Jorgensen 1983；AmberTools 无现成 '
                       'leaprc，getlmp 内置；离子参数沿用 JC_TIP3P 近似）'),
    'tip4pew': dict(leaprc='leaprc.water.tip4pew', box_off='tip4pewbox.off',
                    resname='WAT', nsite=4, supported=True, has_ep=True,
                    angle_k=100.0, angle_theta=104.52, om_dist=0.1250,
                    desc='TIP4P-Ew（4-site，Ewald 优化，Horn 2004）'),
    'tip4pd': dict(leaprc='leaprc.water.tip4pd', box_off='tip4pd.off',
                   resname='WAT', nsite=4, supported=True, has_ep=True,
                   angle_k=100.0, angle_theta=104.52, om_dist=0.1546,
                   desc='TIP4P-D（4-site，色散校正，Piana 2015；'
                        '离子参数 CHARMM22）'),
    'opc': dict(leaprc='leaprc.water.opc', box_off='opcbox.off',
                resname='WAT', nsite=4, supported=True, has_ep=True,
                angle_k=100.0, angle_theta=103.60, om_dist=0.1594,
                desc='OPC（4-site，最优 4 点水，Izadi 2014）'),
}


def water_model(name: str) -> dict:
    """查水模型表（不存在报错，附可用清单）。"""
    m = WATER_MODELS.get(name)
    if m is None:
        raise ValueError(
            f'water.model={name!r} 不支持。可用: '
            f'{sorted(WATER_MODELS)}')
    return m


# ---------------------------------------------------------------- 水坐标
def _parse_off_first_residue(off_path: str,
                             nsite: int = 3) -> list[tuple[str, str, float, float, float]]:
    """解析 off 盒文件第一残基 → [(原子名, 元素, x, y, z), ...]。

    纯文本解析：atoms 表取原子名（第 1 列 "name"），positions 表前
    N 行取坐标（第一残基 = 前 N 行，N = 残基原子数）。
    与 tleap 加载的库严格一致（坐标/原子名/残基模板同一来源）。
    nsite: 取前几个原子（3-site 水取 3；4-site 含 EPW 取 4）。
    """
    txt = open(off_path, encoding='utf-8').read()
    m_atoms = re.search(r'!entry\.\S+\.unit\.atoms table.*?\n(.*?)(?=\n!entry\.|\Z)',
                        txt, re.S)
    if not m_atoms:
        raise RuntimeError(f'off 文件无 atoms 表: {off_path}')
    atom_lines = [l for l in m_atoms.group(1).splitlines() if l.strip()]
    m_pos = re.search(r'!entry\.\S+\.unit\.positions table.*?\n(.*?)(?=\n!entry\.|\Z)',
                      txt, re.S)
    if not m_pos:
        raise RuntimeError(f'off 文件无 positions 表: {off_path}')
    pos_lines = [l for l in m_pos.group(1).splitlines() if l.strip()]

    res_atoms: list[tuple[str, str, float, float, float]] = []
    for i, al in enumerate(atom_lines):
        # 原子行: "name" "type" typex resx flags seq elmnt chg ...
        m = re.match(r'\s*"([^"]+)"\s+"([^"]+)"', al)
        if not m:
            continue
        name = m.group(1)
        elem = _element_from_name(name)
        if i >= len(pos_lines):
            raise RuntimeError(f'off positions 行数不足: {off_path}')
        x, y, z = (float(v) for v in pos_lines[i].split()[:3])
        res_atoms.append((name, elem, x, y, z))
        if len(res_atoms) == nsite:
            break   # 只取第一残基前 nsite 个原子（水 3；4-site 含 EPW 取 4）
    if len(res_atoms) < nsite:
        raise RuntimeError(f'off 第一残基原子数不足 {nsite}: {off_path}')
    return res_atoms


def _element_from_name(name: str) -> str:
    """从原子名推元素（O/H1/H2 → O/H；Na+/Cl- → Na/Cl）。"""
    s = name.rstrip('+-0123456789')
    if s and s[0].isalpha():
        if len(s) >= 2 and s[1].islower():
            return s[:2].capitalize()
        return s[0].upper()
    return s or '?'

# ---------------------------------------------------------------- 离子表
def scan_ion_table() -> dict[str, dict]:
    """扫描 atomic_ions.lib → 单原子离子类型表。

    key = 类型名（如 'Na+'），value = {resname, atomname, elem, charge}。
    - 块结构 `!entry.<块名>.unit.atoms table ...`（块名 = 残基名；
      部分带电荷块另有 `.unit.name` 段，直接用块名即可）；
    - 只保留单原子离子（多原子离子块 atoms 表多行，按需求排除）；
    - 元素从 atoms 表 elmnt 字段（原子序数）取，不靠名字猜
      （旧式全大写残基名 AG/BR/CA 与双字母元素无法用名字可靠区分）；
    - 同一类型可能有多个残基名（如 Ag+ 有 AG / Ag 两块；Na+ 有 Na+ / NA），
      resname 优先取带电荷的块名（Na+ 优于 NA，可读性好）；
    - 原子名取该残基块 atoms 表的原子名（PDB 必须与 lib 模板一致才能被
      tleap loadpdb 匹配）。
    """
    path = os.path.join(LIB_DIR, 'atomic_ions.lib')
    if not os.path.exists(path):
        raise RuntimeError(f'找不到 atomic_ions.lib: {path}（AmberTools 未安装？）')
    txt = open(path, encoding='utf-8').read()
    blocks = re.split(r'(?=!entry\.)', txt)
    types: dict[str, dict] = {}
    for b in blocks:
        m = re.search(
            r'!entry\.([^\s.]+)\.unit\.atoms table.*?\n'
            r'\s*"([^"]+)"\s+"([^"]+)"\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+([\d.+-]+)',
            b, re.S)
        if not m:
            continue
        rname, aname, atype = m.group(1), m.group(2), m.group(3)
        elmnt = int(m.group(4))
        chg = float(m.group(5))
        elem = _element_from_atomic_number(elmnt)
        # 原子数：atoms 表数据行数（atomspertinfo 等行格式不同不会误匹配）
        n_atom_rows = len(re.findall(
            r'"([^"]+)"\s+"([^"]+)"\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+([\d.+-]+)', b))
        if n_atom_rows != 1:
            continue   # 多原子离子排除（H3O+ / NH4+ 等）
        cur = types.get(atype)
        if cur is None:
            types[atype] = dict(resname=rname, atomname=aname,
                                elem=elem, charge=chg)
        else:
            # 同一类型多个残基名：优先带电荷名字（Na+ 优于 NA）
            if _has_charge_mark(rname) and not _has_charge_mark(cur['resname']):
                cur['resname'] = rname
                cur['atomname'] = aname
                cur['charge'] = chg
    return types


def _element_from_atomic_number(z: int) -> str:
    """原子序数 → 元素符号（RDKit 周期表，覆盖全部元素）。"""
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    try:
        return pt.GetElementSymbol(z)
    except Exception:
        return '?'


def _has_charge_mark(name: str) -> bool:
    return '+' in name or '-' in name


def ion_type(name: str) -> dict:
    """查离子类型表（不存在报错，附可用清单）。"""
    table = scan_ion_table()
    t = table.get(name)
    if t is None:
        avail = ', '.join(sorted(table))
        raise ValueError(
            f'ions 中 name={name!r} 不在 AmberTools 单原子离子库中。'
            f'可用类型（atomic_ions.lib 单原子，{len(table)} 个）: {avail}')
    return t

# ---------------------------------------------------------------- 模板 xyz
def water_template_xyz(model: str) -> list[tuple[str, str, float, float, float]]:
    """水模型单分子模板 → [(原子名, 元素, x, y, z), ...]（坐标来自 AmberTools off）。

    4-site 模型返回 4 原子（O/H1/H2/EPW）：EPW 是虚拟位点，坐标保留在模板中
    供 packmol 刚体放置；导出层会丢弃 EPW 并把电荷合并到 O（见 prmtop_to_lammps）。
    """
    m = water_model(model)
    off = os.path.join(LIB_DIR, m['box_off'])
    if not os.path.exists(off):
        raise RuntimeError(f'找不到水模型盒文件: {off}（AmberTools 数据缺失）')
    return _parse_off_first_residue(off, nsite=m['nsite'])


def water_om_distance(model: str) -> float:
    """4-site 水 O→M(EPW) 距离（Å），从 off 模板坐标实测。

    与 tleap 加载的库几何严格一致；LAMMPS pair_style lj/cut/tip4p/long 的
    qdist 参数用此值。3-site 模型返回 None。
    """
    rows = water_template_xyz(model)
    names = [r[0] for r in rows]
    if 'EPW' not in names:
        return None
    io, ie = names.index('O'), names.index('EPW')
    o, e = rows[io], rows[ie]
    return sum((o[k] - e[k]) ** 2 for k in range(2, 5)) ** 0.5


def ion_template_xyz(ion: str) -> list[tuple[str, str, float, float, float]]:
    """离子单分子模板 → [(原子名, 元素, 0, 0, 0)]（坐标由 packmol 放置）。"""
    t = ion_type(ion)
    return [(t['atomname'], t['elem'], 0.0, 0.0, 0.0)]


def write_template_xyz(rows: list[tuple[str, str, float, float, float]],
                       out_xyz: str) -> int:
    """模板原子行 → packmol 输入 xyz（元素 + 坐标）。返回原子数。"""
    lines = [str(len(rows)), 'from getlmp (template)']
    for name, elem, x, y, z in rows:
        lines.append(f'{elem:2s} {x:12.5f} {y:12.5f} {z:12.5f}')
    with open(out_xyz, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return len(rows)


def element_mass(elem: str) -> float:
    """元素相对原子质量（RDKit 周期表，覆盖全部元素；密度盒估算用）。

    'E' = 4-site 水虚拟位点 EPW 在 packmol xyz 里的元素列标识（无质量）→ 0。
    """
    if elem == 'E':
        return 0.0
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    try:
        return pt.GetAtomicWeight(elem)
    except Exception:
        raise RuntimeError(f'未知元素 {elem!r}（RDKit 周期表无法识别，密度盒无法估算）')
