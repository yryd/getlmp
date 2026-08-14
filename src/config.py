#!/usr/bin/env python
"""配置层：解析 input.yaml → cfg dict，并做合法性校验。

python main.py input.yaml 的唯一入口配置。支持字段见 docs/。
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import yaml

SUPPORTED_CHARGE = {'bcc', 'abcg2', 'resp', 'resp2', 'none'}
SUPPORTED_FF = {'gaff2', 'gaff', 'reaxff'}
SUPPORTED_QM_ENGINE = {'gaussian', 'quick'}
# ReaxFF 原子类型默认元素顺序（须与 ffield.reax 的元素顺序一致；可经 yaml 覆盖）
DEFAULT_REAX_ELEMENTS = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']


@dataclass
class MoleculeCfg:
    smiles: str = ''         # SMILES（与 xyz 二选一）
    name: str = ''
    count: int = 1
    resname: str = ''        # 残基名（antechamber -rn），默认取 name 前 3 字符大写
    xyz: str = ''            # xyz 输入路径（与 smiles 二选一；绝对路径或相对 yaml）
                             # 原子坐标原样使用（无加 H/无优化），键由 RDKit 推断


@dataclass
class WaterCfg:
    """水溶剂（内置模板，AmberTools 原生力场）。

    yaml 写法：
      water:
        model: tip3p        # 一期: tip3p / spce / opc3（二期: opc/tip4pew/tip4pd 预留）
        count: 3000         # 水分子数（数量，不做浓度）
    出现即强制多分子 packmol 装盒。
    """
    model: str = ''
    count: int = 0


@dataclass
class IonCfg:
    """离子（AmberTools atomic_ions.lib 单原子库，LJ 随 leaprc 加载）。

    yaml 写法：
      ions:
        - name: Na+         # atomic_ions.lib 类型名（Na+ Cl- K+ Ca2+ ... 单原子）
          count: 20
    多原子离子（H3O+/NH4+ 等）不支持。
    """
    name: str = ''
    count: int = 0


@dataclass
class PackmolCfg:
    enabled: bool = False    # count>1 时自动启用
    preset: str = 'bulk'     # bulk / slab / interface（阶段 2 先实现 bulk）
    box: list = None         # [xlo,ylo,zlo,xhi,yhi,zhi]（与 density 二选一；density 优先）
    seed: int = 2026
    tolerance: float = 2.0   # Å
    nloop0: int = 1000       # packmol 初始随机放置尝试次数（默认 20 对高密度大体系常不够，报错提示加大）
    density: float = 0.0     # 目标密度 g/cm³（>0 时按 总质量/密度 自动算立方盒，忽略 box）
    inp_file: str = ''       # 自定义 packmol inp（绝对路径，load_config 已解析）。
                             # 非空时跳过自动生成 write_inp，直接用该文件跑 packmol；
                             # structure 顺序须与 molecules 一致（可改 number 与约束行）。
                             # 与 density 互斥（自动盒子只对自动生成的 inp 生效）。


@dataclass
class EspCfg:
    """ESP 可视化导出（单分子 RESP/RESP2 时默认生成 iso 产物，pt 可选）。

    产物位于 _others/electrostatic_potential/ 下：
      iso/{density.cub, esp.cub}   iso 法（VMD: 密度等值面 + ESP 着色）
      pt/{mol.pdb, vtx.pdb}        pt 法（分子结构 + ESP 曲面顶点，Beta 着色）
    统一由 Multiwfn 导出（gaussian 用 .fch，quick 用 .molden）。

    spacing/timeout 支持 'auto'（按原子数自动分档，见 multiwfn.auto_esp_params）
    或显式指定。2026-08-14：Multiwfn iso 法默认表面格点间距 0.25 Å 对
    ≥1000 基函数大分子点数爆炸超时，故引入自动分档。
    """
    enabled: bool = True     # 仅单分子 RESP/RESP2 生效
    pt: bool = False         # pt 法产物（mol.pdb/vtx.pdb），默认关（基本不用；需要时开）
    spacing: object = 'auto'     # 'auto'=按原子数分档（≤20:0.25 / 21-40:0.3 / >40:0.4）| 0.15~0.8 显式
    timeout: object = 'auto'     # 'auto'=按原子数分档（≤20:600 / 21-40:1800 / >40:3600）| 秒数；0=不限
    # 旧字段 spacing/buffer（废弃点电荷 cube 参数）已移除；yaml 里遗留的 buffer 键被忽略


@dataclass
class QmCfg:
    """量子化学引擎配置（RESP/RESP2 电荷拟合的波函数来源）。

    engine=gaussian（默认，新主路径）：G16 单点 + Multiwfn RESP 拟合；
    engine=quick（旧路径保留）：QUICK HF/6-31G* + resp 两阶段拟合。
    """
    engine: str = 'gaussian'   # gaussian / quick
    g16root: str = ''          # 显式 G16 根目录覆盖（一般留空，PATH 已含 g16 即可）
    method: str = 'b3lyp'      # 531 推荐 B3LYP-D3(BJ)
    basis: str = 'def2TZVP'    # G16 命名（def2-TZVP 连字符写法会自动转）
    opt: bool = False          # 是否先几何优化（默认单点）
    solvent: str = ''          # 空=气相；'water'/'ethanol'（RESP2 需要）
    resp2: bool = False        # RESP2（需 solvent 非空）
    delta: float = 0.5         # RESP2 δ
    multiwfn_path: str = ''    # 空=自动探测


@dataclass
class Config:
    name: str = 'system'
    forcefield: str = 'gaff2'
    charge_method: str = 'bcc'
    net_charge: int = 0
    molecules: list = field(default_factory=list)   # list[MoleculeCfg]
    water: WaterCfg | None = None                    # 水溶剂（None=无水；出现即多分子）
    ions: list = field(default_factory=list)         # list[IonCfg]
    packmol: PackmolCfg = field(default_factory=PackmolCfg)
    esp: EspCfg = field(default_factory=EspCfg)
    qm: QmCfg = field(default_factory=QmCfg)
    output: str = 'data.lmp'
    workdir: str = 'work'
    seed: int = 2026
    organize_output: bool = True   # 跑完只留 lmp+mol2，其余进 workdir/_others/
    organize_backup: bool = False  # _others/ 已存在同名时：false=覆盖只留最新，true=追加 HHMMSS 时间戳保留多份
    reuse_molecule: bool = False   # 分子层复用：workdir 已有同配置指纹的 mol2/frcmod/波函数时跳过
                                   # 重算（默认关）。改 SMILES/电荷方法/力场/qm 后指纹不匹配会自动重算。
    buffer: float = 3.8        # 单分子盒 padding（Å，仅 count=1 时使用）
    reax_elements: list = field(default_factory=lambda: list(DEFAULT_REAX_ELEMENTS))
    reax_atom_style: str = 'charge'   # charge(6 列,无 mol-id) / full(7 列,带 mol-id)

    @property
    def base_dir(self) -> str:
        """配置文件所在目录（相对路径以此为准）。"""
        return os.path.dirname(self._source_path) if self._source_path else '.'

    _source_path: str = ''


def load_config(path: str) -> Config:
    """读取 input.yaml → Config（校验字段与取值）。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f'找不到配置文件: {path}')
    with open(path, encoding='utf-8') as f:
        raw = yaml.safe_load(f) or {}

    cfg = Config()
    cfg._source_path = os.path.abspath(path)

    cfg.name = str(raw.get('name', cfg.name))
    ff = raw.get('forcefield', cfg.forcefield)
    if ff not in SUPPORTED_FF:
        raise ValueError(f'不支持的 forcefield={ff!r}（支持: {sorted(SUPPORTED_FF)}）')
    cfg.forcefield = ff

    # qm 段（RESP/RESP2 的量子化学引擎配置；bcc/abcg2 不依赖 qm）
    qm = raw.get('qm') or {}
    if not isinstance(qm, dict):
        raise ValueError('qm 需要字典配置，如 {engine: gaussian, method: b3lyp, basis: def2-TZVP}')
    qe = str(qm.get('engine', cfg.qm.engine))
    if qe not in SUPPORTED_QM_ENGINE:
        raise ValueError(f'不支持的 qm.engine={qe!r}（支持: {sorted(SUPPORTED_QM_ENGINE)}）')
    cfg.qm.engine = qe
    cfg.qm.g16root = str(qm.get('g16root', cfg.qm.g16root))
    cfg.qm.method = str(qm.get('method', cfg.qm.method))
    cfg.qm.basis = str(qm.get('basis', cfg.qm.basis))
    cfg.qm.opt = bool(qm.get('opt', cfg.qm.opt))
    cfg.qm.solvent = str(qm.get('solvent', cfg.qm.solvent))
    cfg.qm.resp2 = bool(qm.get('resp2', cfg.qm.resp2))
    cfg.qm.delta = float(qm.get('delta', cfg.qm.delta))
    cfg.qm.multiwfn_path = str(qm.get('multiwfn_path', cfg.qm.multiwfn_path))
    if not (0.0 <= cfg.qm.delta <= 1.0):
        raise ValueError(f'qm.delta 需在 [0,1] 内，当前 {cfg.qm.delta!r}')

    cm = raw.get('charge_method', cfg.charge_method)
    if cm not in SUPPORTED_CHARGE:
        raise ValueError(
            f'不支持的 charge_method={cm!r}（支持: {sorted(SUPPORTED_CHARGE)}）')
    cfg.charge_method = cm
    if cm in ('resp', 'resp2'):
        if qe == 'quick' and cm == 'resp2':
            raise ValueError('charge_method=resp2 需要 qm.engine=gaussian'
                             '（QUICK 无 resp2 方法）')
        if cm == 'resp2' and not cfg.qm.resp2:
            raise ValueError('charge_method=resp2 需要 qm.resp2: true')
        if cm == 'resp2' and not cfg.qm.solvent:
            raise ValueError('charge_method=resp2 需要 qm.solvent 非空'
                             '（RESP2 = 气相 + 溶剂 PCM 单点混合）')
    if ff == 'reaxff' and cm not in ('none',):
        raise ValueError(
            f'forcefield=reaxff 时 charge_method 应为 none（QEq 在模拟中计算），'
            f'当前 {cm!r}')

    # ReaxFF 元素顺序（决定 data 类型号，须与 ffield.reax 一致）
    rx = raw.get('reax_elements')
    if rx:
        if not isinstance(rx, list) or not all(isinstance(e, str) for e in rx):
            raise ValueError('reax_elements 需要字符串列表，如 [C, H, O]')
        seen = []
        for e in rx:
            if e not in seen:
                seen.append(e)
        cfg.reax_elements = seen

    if ff == 'reaxff':
        ras = raw.get('reax_atom_style', cfg.reax_atom_style)
        if ras not in ('charge', 'full'):
            raise ValueError(f'reax_atom_style 只支持 charge/full，当前 {ras!r}')
        cfg.reax_atom_style = ras

    cfg.net_charge = int(raw.get('net_charge', cfg.net_charge))
    cfg.output = str(raw.get('output', cfg.output))
    cfg.workdir = str(raw.get('workdir', cfg.workdir))
    cfg.seed = int(raw.get('seed', cfg.seed))
    cfg.organize_output = bool(raw.get('organize_output', cfg.organize_output))
    cfg.organize_backup = bool(raw.get('organize_backup', cfg.organize_backup))
    cfg.reuse_molecule = bool(raw.get('reuse_molecule', cfg.reuse_molecule))
    cfg.buffer = float(raw.get('buffer', cfg.buffer))

    mols = raw.get('molecules')
    if not mols:
        # 纯溶剂/离子体系允许无 molecules（water 或 ions 存在时）
        if not raw.get('water') and not raw.get('ions'):
            raise ValueError('配置缺少 molecules 段（或 water/ions 溶剂段）')
        mols = []
    multi = False
    for m in mols:
        smiles = m.get('smiles', '')
        xyz = str(m.get('xyz', '') or '').strip()
        if not smiles and not xyz:
            raise ValueError('molecules 条目需要 smiles 或 xyz（二选一）')
        if smiles and xyz:
            raise ValueError('molecules 条目不能同时给 smiles 和 xyz（二选一）')
        name = str(m.get('name', f'mol{len(cfg.molecules) + 1}'))
        mc = MoleculeCfg(
            smiles=str(smiles),
            name=name,
            count=int(m.get('count', 1)),
            resname=str(m.get('resname', name[:3].upper())),
        )
        if xyz:
            p = xyz if os.path.isabs(xyz) else os.path.join(cfg.base_dir, xyz)
            if not os.path.exists(p):
                raise ValueError(f'molecules 的 xyz 文件不存在: {p}')
            mc.xyz = os.path.abspath(p)
        if mc.count > 1:
            multi = True
        cfg.molecules.append(mc)

    # water 段（内置水模板；出现即强制多分子装盒）
    wt = raw.get('water')
    if wt is not None:
        if not isinstance(wt, dict):
            raise ValueError('water 需要字典配置，如 {model: tip3p, count: 3000}')
        from solvent_templates import water_model   # 查表（含二期拦截）
        model = str(wt.get('model', '')).strip()
        if not model:
            raise ValueError('water.model 不能为空（一期: tip3p / spce / opc3）')
        wm = water_model(model)   # 不存在/二期预留 → 报错
        count = int(wt.get('count', 0))
        if count <= 0:
            raise ValueError(f'water.count 需为正整数（分子数），当前 {count!r}')
        cfg.water = WaterCfg(model=model, count=count)
        multi = True
        print(f'  [config] 水溶剂: {model} × {count}（模板 {wm["box_off"]}）')

    # ions 段（atomic_ions.lib 单原子；出现即强制多分子装盒）
    ions = raw.get('ions')
    if ions:
        if not isinstance(ions, list):
            raise ValueError('ions 需要列表，如 [{name: Na+, count: 20}]')
        from solvent_templates import ion_type   # 查表（单原子库）
        for it in ions:
            if not isinstance(it, dict):
                raise ValueError(f'ions 条目需要字典 {{name, count}}，收到 {it!r}')
            iname = str(it.get('name', '')).strip()
            if not iname:
                raise ValueError('ions 条目缺少 name（如 Na+ / Cl-）')
            t = ion_type(iname)   # 不在库 → 报错（附可用清单）
            icount = int(it.get('count', 0))
            if icount <= 0:
                raise ValueError(f'ions 中 {iname} 的 count 需为正整数，当前 {icount!r}')
            cfg.ions.append(IonCfg(name=iname, count=icount))
            multi = True
            print(f'  [config] 离子: {iname} × {icount}'
                  f'（atomic_ions.lib, 电荷 {t["charge"]:+.0f}）')

    if cfg.forcefield == 'reaxff' and (cfg.water is not None or cfg.ions):
        raise ValueError('forcefield=reaxff 暂不支持内置水/离子模板'
                         '（ReaxFF 的水/离子请走 molecules 输入；'
                         '内置模板仅 gaff/gaff2 力场）')

    # packmol 段（多分子自动启用 bulk）
    pm = raw.get('packmol') or {}
    cfg.packmol.preset = str(pm.get('preset', 'bulk'))
    box = pm.get('box')
    if box:
        if len(box) != 6:
            raise ValueError('packmol.box 需要 6 个数 [xlo,ylo,zlo,xhi,yhi,zhi]')
        cfg.packmol.box = [float(v) for v in box]
    cfg.packmol.seed = int(pm.get('seed', cfg.seed))
    cfg.packmol.tolerance = float(pm.get('tolerance', 2.0))
    cfg.packmol.nloop0 = int(pm.get('nloop0', 1000))
    cfg.packmol.density = float(pm.get('density', 0.0))
    cfg.packmol.enabled = multi or bool(pm.get('enabled', False))
    if cfg.packmol.enabled and cfg.packmol.box is None and cfg.packmol.density <= 0:
        raise ValueError('多分子体系需要 packmol.box（如 [0,0,0,60,60,60]）'
                         '或 packmol.density（g/cm³，自动按总质量算立方盒）')
    if cfg.packmol.density <= 0:
        cfg.packmol.density = 0.0
    if cfg.packmol.enabled and cfg.packmol.preset != 'bulk':
        raise ValueError(f'阶段 2 仅支持 packmol preset=bulk，当前 {cfg.packmol.preset!r}；'
                         f'slab/interface 见规划文档（后续阶段）')
    inp_file = str(pm.get('inp_file', '')).strip()
    if inp_file:
        if not cfg.packmol.enabled:
            raise ValueError('packmol.inp_file 仅多分子体系（count>1 或 packmol.enabled）可用')
        if cfg.packmol.density > 0:
            raise ValueError('packmol.inp_file 与 packmol.density 互斥'
                             '（自动盒子只对自动生成的 inp 生效；自定义 inp 里可写 density 关键字，'
                             '但盒子由 packmol 内部计算，data.lmp 盒边界无法对齐，故不支持）')
        p = inp_file if os.path.isabs(inp_file) else os.path.join(cfg.base_dir, inp_file)
        if not os.path.exists(p):
            raise ValueError(f'packmol.inp_file 不存在: {p}')
        cfg.packmol.inp_file = os.path.abspath(p)

    # ESP 可视化导出（默认开启；仅单分子 RESP/RESP2 生效）
    esp = raw.get('esp') or {}
    if not isinstance(esp, dict):
        raise ValueError('esp 需要字典配置，如 {enabled: true, spacing: auto, timeout: auto}')
    cfg.esp.enabled = bool(esp.get('enabled', True))
    cfg.esp.pt = bool(esp.get('pt', False))
    # spacing: 'auto'（按原子数分档）或 0.15~0.8 Å 显式指定
    s = esp.get('spacing', 'auto')
    if isinstance(s, str):
        if s.lower() != 'auto':
            raise ValueError(f'esp.spacing 只支持 auto 或数值（如 0.3），收到 {s!r}')
        cfg.esp.spacing = 'auto'
    else:
        s = float(s)
        if not (0.15 <= s <= 0.8):
            raise ValueError(f'esp.spacing 需在 0.15~0.8 Å（Multiwfn 手册建议范围），收到 {s}')
        cfg.esp.spacing = s
    # timeout: 'auto'（按原子数分档）或秒数；0=不限（subprocess timeout=None）
    t = esp.get('timeout', 'auto')
    if isinstance(t, str):
        if t.lower() != 'auto':
            raise ValueError(f'esp.timeout 只支持 auto、秒数或 0（不限），收到 {t!r}')
        cfg.esp.timeout = 'auto'
    else:
        t = int(t)
        if t < 0:
            raise ValueError(f'esp.timeout 不能为负数: {t}')
        cfg.esp.timeout = t
    return cfg
