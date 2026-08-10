#!/usr/bin/env python
"""配置层：解析 input.yaml → cfg dict，并做合法性校验。

python main.py input.yaml 的唯一入口配置。支持字段见 docs/。
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import yaml

SUPPORTED_CHARGE = {'bcc', 'abcg2', 'resp', 'none'}   # resp2 不支持（见 docs/dev_notes.md）
SUPPORTED_FF = {'gaff2', 'gaff', 'reaxff'}
# ReaxFF 原子类型默认元素顺序（须与 ffield.reax 的元素顺序一致；可经 yaml 覆盖）
DEFAULT_REAX_ELEMENTS = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']


@dataclass
class MoleculeCfg:
    smiles: str
    name: str
    count: int = 1
    resname: str = ''        # 残基名（antechamber -rn），默认取 name 前 3 字符大写


@dataclass
class PackmolCfg:
    enabled: bool = False    # count>1 时自动启用
    preset: str = 'bulk'     # bulk / slab / interface（阶段 2 先实现 bulk）
    box: list = None         # [xlo,ylo,zlo,xhi,yhi,zhi]
    seed: int = 2026
    tolerance: float = 2.0   # Å


@dataclass
class EspCfg:
    """ESP 可视化导出（单分子 RESP 时默认不生成，需要时 esp.enabled: true 开启）。"""
    enabled: bool = False    # 默认关闭；仅单分子 RESP 生效
    spacing: float = 0.3     # 网格间距 Å（0.2 更细、0.5 更快）
    buffer: float = 1.5      # 分子外扩 Å（网格覆盖范围）


@dataclass
class Config:
    name: str = 'system'
    forcefield: str = 'gaff2'
    charge_method: str = 'bcc'
    net_charge: int = 0
    molecules: list = field(default_factory=list)   # list[MoleculeCfg]
    packmol: PackmolCfg = field(default_factory=PackmolCfg)
    esp: EspCfg = field(default_factory=EspCfg)
    output: str = 'data.lmp'
    workdir: str = 'work'
    seed: int = 2026
    organize_output: bool = True   # 跑完只留 lmp+mol2，其余进 workdir/_others/
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

    cm = raw.get('charge_method', cfg.charge_method)
    if cm not in SUPPORTED_CHARGE:
        raise ValueError(
            f'不支持的 charge_method={cm!r}（支持: {sorted(SUPPORTED_CHARGE)}；'
            f'resp2 不支持：antechamber 无该电荷方法（见 docs/dev_notes.md）')
    cfg.charge_method = cm
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
    cfg.buffer = float(raw.get('buffer', cfg.buffer))

    mols = raw.get('molecules')
    if not mols:
        raise ValueError('配置缺少 molecules 段')
    multi = False
    for m in mols:
        smiles = m.get('smiles')
        if not smiles:
            raise ValueError('molecules 条目缺少 smiles')
        name = str(m.get('name', f'mol{len(cfg.molecules) + 1}'))
        mc = MoleculeCfg(
            smiles=str(smiles),
            name=name,
            count=int(m.get('count', 1)),
            resname=str(m.get('resname', name[:3].upper())),
        )
        if mc.count > 1:
            multi = True
        cfg.molecules.append(mc)

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
    cfg.packmol.enabled = multi or bool(pm.get('enabled', False))
    if cfg.packmol.enabled and cfg.packmol.box is None:
        raise ValueError('多分子体系需要 packmol.box（如 [0,0,0,60,60,60]）')
    if cfg.packmol.enabled and cfg.packmol.preset != 'bulk':
        raise ValueError(f'阶段 2 仅支持 packmol preset=bulk，当前 {cfg.packmol.preset!r}；'
                         f'slab/interface 见规划文档（后续阶段）')

    # ESP 可视化导出（默认关闭；仅单分子 RESP 生效，需要时 esp.enabled: true 开启）
    esp = raw.get('esp') or {}
    if not isinstance(esp, dict):
        raise ValueError('esp 需要字典配置，如 {enabled: true, spacing: 0.3, buffer: 1.5}')
    cfg.esp.enabled = bool(esp.get('enabled', False))
    cfg.esp.spacing = float(esp.get('spacing', 0.3))
    cfg.esp.buffer = float(esp.get('buffer', 1.5))
    if cfg.esp.spacing <= 0 or cfg.esp.buffer <= 0:
        raise ValueError('esp.spacing / esp.buffer 需为正数')
    return cfg
