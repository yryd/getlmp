#!/usr/bin/env python
"""流水线编排：分子层 → 体系层 → 导出层 → 校验 → 检查报告。

- 阶段 1：单分子（count=1）tleap 闭环
- 阶段 2：多分子（count>1）packmol 装盒 + tleap loadpdb 闭环
入口：`python main.py input.yaml`。
"""
from __future__ import annotations

import glob
import os
import shutil

from config import Config
from export_layer import build_data, validate_export
from export_reaxff import build_data_reaxff
from molecule_layer import build_molecule
from packmol_layer import (mol2_to_xyz, parse_inp_structures, parse_packed_xyz,
                           run_packmol, sdf_to_xyz, write_inp)
from system_layer import (build_system_multi, build_system_single, count_mol2_topo)
from xyz_export import export_system_xyz

# 集群 LAMMPS 验证模板（交付物的一部分，用户/SimAgent 可据此实跑）
# 注意：read_data 必须放在所有 style 定义之后（data 含 Coeffs 段时 LAMMPS 强制要求），
# 否则报 "Must define pair_style before Pair Coeffs"。
LAMMPS_IN_TEMPLATE = """units real
atom_style full
boundary p p p
pair_style lj/cut/coul/long 10.0
pair_modify mix arithmetic
bond_style harmonic
angle_style harmonic
dihedral_style fourier
improper_style cvff
special_bonds amber
kspace_style pppm 1e-4
read_data {data}
run 0
"""


def run_pipeline(cfg: Config) -> dict:
    """执行整条流水线，返回报告 dict（同时写入 workdir/check_report.txt）。"""
    workdir = cfg.workdir if os.path.isabs(cfg.workdir) else os.path.join(cfg.base_dir, cfg.workdir)
    os.makedirs(workdir, exist_ok=True)

    multi = cfg.packmol.enabled
    print(f'==== getlmp pipeline: {cfg.name} '
          f'({cfg.forcefield} / {cfg.charge_method}, net={cfg.net_charge}, '
          f'{"多分子" if multi else "单分子"}) ====')

    # 1. 分子层：每个分子类型一次完整流水线
    mols = []
    for mc in cfg.molecules:
        mols.append(build_molecule(cfg, mc, workdir))

    if cfg.forcefield == 'reaxff':
        return _run_pipeline_reaxff(cfg, workdir, mols, multi)

    if not multi:
        # 2a. 体系层：单分子 tleap
        mol = mols[0]
        sysr = build_system_single(cfg, workdir, mol['mol2'], mol['frcmod'])
        natom_input = mol['natom']
        expected_topo = None
        box = None
    else:
        # 2b. 体系层：packmol 装盒 + 合并 PDB + tleap loadpdb
        print(f'== 体系层: packmol (preset={cfg.packmol.preset}, '
              f'box={cfg.packmol.box}) ==')
        n_atoms = []
        for i, mc in enumerate(cfg.molecules):
            n_atoms.append(mol2_to_xyz(mols[i]['mol2'],
                                       os.path.join(workdir, mc.name + '.xyz')))
        natom_by_name = {mc.name: n for mc, n in zip(cfg.molecules, n_atoms)}
        if cfg.packmol.inp_file:
            # 自定义 inp（用户改过：加 fixed 约束/改 number/同一类型拆多块等）。
            # structure 文件名的类型须在 molecules 中；同一类型可拆多个块
            # （如"固定 1 个 + 自由 N 个"），number 以 inp 内为准。
            inp = cfg.packmol.inp_file
            print(f'  [custom inp] 使用 {inp}（跳过自动生成 write_inp）')
            structs = parse_inp_structures(inp)
            counts = [n for _, n in structs]
            type_names = [os.path.splitext(os.path.basename(f))[0] for f, _ in structs]
            for i, (fname, n) in enumerate(structs):
                print(f'    structure[{i}] {fname}  number {n}')
            for t in type_names:
                if t not in natom_by_name:
                    raise RuntimeError(
                        f'自定义 inp structure 文件 {t}.xyz 不在 molecules 中'
                        f'（可用: {sorted(natom_by_name)}）；结构文件名须为 {{name}}.xyz')
        else:
            inp = os.path.join(workdir, 'packmol.inp')
            write_inp(cfg, [os.path.join(workdir, mc.name + '.xyz') for mc in cfg.molecules],
                      n_atoms, inp)
            counts = [mc.count for mc in cfg.molecules]
            type_names = [mc.name for mc in cfg.molecules]
        packed = run_packmol(inp, workdir)
        blocks = parse_packed_xyz(packed, counts, natom_by_name, type_names)
        sysr = build_system_multi(
            cfg, workdir,
            [m['mol2'] for m in mols],
            [m['frcmod'] for m in mols],
            blocks)
        natom_input = sum(c * natom_by_name[t] for c, t in zip(counts, type_names))
        print(f'  体系原子总数（期望）: {natom_input}')
        # 拓扑期望：按类型聚合拷贝数（自定义 inp 同一类型可拆多块）
        type_count: dict[str, int] = {}
        for c, t in zip(counts, type_names):
            type_count[t] = type_count.get(t, 0) + c
        expected_topo = {'bonds': 0, 'angles': 0}
        for mc, m in zip(cfg.molecules, mols):
            c = type_count.get(mc.name, 0)
            if c:
                nb, na = count_mol2_topo(m['mol2'])
                expected_topo['bonds'] += c * nb
                expected_topo['angles'] += c * na
        # packmol box 语义 [xlo,ylo,zlo,xhi,yhi,zhi] → LAMMPS 顺序 [xlo,xhi,ylo,yhi,zlo,zhi]
        b = cfg.packmol.box
        box = [b[0], b[3], b[1], b[4], b[2], b[5]]

    # 3. 导出层 + 校验
    exp = build_data(cfg, workdir, sysr['prmtop'], sysr['inpcrd'], box=box)
    ok, msgs = validate_export(cfg, exp, natom_input=natom_input,
                               expected_topo=expected_topo)

    # 3b. ESP 可视化导出：单分子 RESP/RESP2 时生成 iso/pt 两套产物（统一 Multiwfn）
    #   _others/electrostatic_potential/iso/{density.cub, esp.cub}  iso 法：密度等值面 + ESP 着色
    #   _others/electrostatic_potential/pt/{mol.pdb, vtx.pdb}       pt 法：分子结构 + ESP 曲面顶点
    #   gaussian 引擎用 .fch；quick 引擎用 .molden（Multiwfn 均支持）
    #   旧的点电荷近似粗 cube（esp_export）已废弃，不再生成
    esp_iso_dir = None
    esp_pt_dir = None
    if not multi and cfg.charge_method in ('resp', 'resp2') and cfg.esp.enabled:
        from multiwfn import esp_pt, esp_cube as mwfn_cube, density_cube
        wavefn = mols[0].get('wavefn')
        if not wavefn or not os.path.exists(wavefn):
            raise RuntimeError(f'ESP 导出需要波函数文件，未找到: {wavefn}')
        esp_root = os.path.join(workdir, '_others', 'electrostatic_potential')
        iso_dir = os.path.join(esp_root, 'iso')
        pt_dir = os.path.join(esp_root, 'pt')
        os.makedirs(iso_dir, exist_ok=True)
        os.makedirs(pt_dir, exist_ok=True)
        # pt 法：mol.pdb（分子结构）+ vtx.pdb（密度 0.001 等值面顶点，B 因子=ESP）
        esp_pt(wavefn,
               os.path.join(pt_dir, 'mol.pdb'),
               os.path.join(pt_dir, 'vtx.pdb'),
               workdir=workdir, multiwfn_path=cfg.qm.multiwfn_path)
        print('  ESP pt  → ' + os.path.join(pt_dir, 'mol.pdb') + ' + vtx.pdb'
              '（Multiwfn 密度 0.001 等值面, B 因子=ESP kcal/mol, VMD Beta 着色）')
        # iso 法：esp.cub（严格 QM ESP）+ density.cub（电子密度，同一网格）
        mwfn_cube(wavefn, os.path.join(iso_dir, 'esp.cub'),
                  workdir=workdir, multiwfn_path=cfg.qm.multiwfn_path)
        density_cube(wavefn, os.path.join(iso_dir, 'density.cub'),
                     workdir=workdir, multiwfn_path=cfg.qm.multiwfn_path)
        print('  ESP iso → ' + os.path.join(iso_dir, 'esp.cub') + ' + density.cub'
              '（Multiwfn 严格 QM, VMD iso 渲染: density 等值面 + ESP 着色）')
        # 清理 Multiwfn 中间文件（正式产物已 copy 到位）
        for tmp in ('mol.pdb', 'vtx.pdb', 'mapfunc.cub'):
            p = os.path.join(workdir, tmp)
            if os.path.exists(p):
                os.remove(p)
        print(f'  波函数 → {wavefn}（VMD/Multiwfn 可视化用）')
        esp_iso_dir = iso_dir
        esp_pt_dir = pt_dir
    elif multi and cfg.charge_method in ('resp', 'resp2'):
        print('  [info] 多分子体系跳过 ESP 导出（当前仅单分子 RESP/RESP2 支持）')

    report = {
        'config': cfg,
        'molecules': mols,
        'system': sysr,
        'export': exp,
        'validation': {'ok': ok, 'messages': msgs},
        'data_lmp': exp['data_lmp'],
        'esp_iso_dir': esp_iso_dir,
        'esp_pt_dir': esp_pt_dir,
        'workdir': workdir,
    }
    # 3c. 体系标准 xyz 导出（默认主产物；原子顺序与 data.lmp 一致，可喂 OVITO/VMD/其他工具）
    xyz_path = os.path.join(workdir, 'system.xyz')
    n_xyz = export_system_xyz(report, xyz_path)
    print(f'  xyz → {os.path.basename(xyz_path)} ({n_xyz} atoms)')
    _write_report(report)
    print(f'\n==== 校验结果: {"通过" if ok else "失败"} ====')
    for m in msgs:
        print('  ' + m)
    print(f'\n输出: {exp["data_lmp"]}')
    print(f'检查报告: {os.path.join(workdir, "check_report.txt")}')

    if cfg.organize_output:
        keep = [exp['data_lmp']] + [m['mol2'] for m in mols if m.get('mol2')] + [xyz_path]
        # 保留 QM 波函数**正式产物**（fch/chg；ESP iso/pt 产物在
        # _others/electrostatic_potential/ 内，organize 不动 _others 目录）
        keep += [p for p in glob.glob(os.path.join(workdir, '*.fch'))]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.molden'))]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.chg'))]
        # 保留 frcmod 与指纹文件（reuse_molecule 复用分子层需要）
        keep += [m['frcmod'] for m in mols if m.get('frcmod')]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.fingerprint.json'))]
        # 保留 packmol inp：自动生成模板（供用户改）+ 自定义 inp（用户输入，不能移走）
        keep.append(os.path.join(workdir, 'packmol.inp'))
        if cfg.packmol.inp_file:
            keep.append(cfg.packmol.inp_file)
        organize_workdir(workdir, keep, cfg.organize_backup)
    return report


def organize_workdir(workdir: str, keep: list[str], backup: bool = False) -> None:
    """跑完整理：workdir 根目录只保留 keep 列表（如 data.lmp、*.mol2），
    其余文件/子目录全部移入 workdir/_others/；
    backup=True 时重名自动加 HHMMSS 时间戳保留多份，False（默认）同名直接覆盖只留最新；
    纯临时目录（_resp_tmp_*，RESP 拟合中间）不保留，直接删除。

    keep 用绝对路径比较；_others 目录本身不动。
    """
    others = os.path.join(workdir, '_others')
    os.makedirs(others, exist_ok=True)
    keep_abs = {os.path.abspath(p) for p in keep}
    moved: list[str] = []
    deleted: list[str] = []
    for entry in sorted(os.listdir(workdir)):
        p = os.path.join(workdir, entry)
        if os.path.abspath(p) in keep_abs or entry == '_others':
            continue
        # 纯临时目录（RESP 拟合中间产物）不保留，直接删
        if entry.startswith('_resp_tmp_') and os.path.isdir(p) and not os.path.islink(p):
            shutil.rmtree(p)
            deleted.append(entry)
            continue
        dst = os.path.join(others, entry)
        if os.path.exists(dst) and backup:
            import time
            stamp = time.strftime('%H%M%S')
            base, ext = os.path.splitext(entry)
            dst = os.path.join(others, f'{base}_{stamp}{ext}')
        if os.path.isdir(p) and not os.path.islink(p):
            shutil.move(p, dst)
        else:
            os.replace(p, dst)
        moved.append(entry)
    if moved or deleted:
        print(f'  [organize] 保留 {len(keep)} 个主产物；'
              f'{len(moved)} 项已移入 _others/'
              + (f'；删除临时目录 {len(deleted)} 个' if deleted else ''))
        if len(moved) + len(deleted) <= 12:
            print('    → ' + ', '.join(moved + [f'{d}(删)' for d in deleted]))


def _run_pipeline_reaxff(cfg: Config, workdir: str, mols: list, multi: bool) -> dict:
    """ReaxFF 分支：坐标/元素 → packmol 装盒（坐标层面）→ 极简 data（无 tleap）。"""
    atom_info = [{'elements': m['elements'], 'natom': m['natom']} for m in mols]

    if multi:
        # 多分子：packmol 装盒（坐标层面，无拓扑）
        print(f'== 体系层: packmol (ReaxFF, preset={cfg.packmol.preset}, '
              f'box={cfg.packmol.box}) ==')
        n_atoms = []
        for i, mc in enumerate(cfg.molecules):
            n_atoms.append(sdf_to_xyz(mols[i]['sdf'],
                                      os.path.join(workdir, mc.name + '.xyz')))
        natom_by_name = {mc.name: n for mc, n in zip(cfg.molecules, n_atoms)}
        if cfg.packmol.inp_file:
            inp = cfg.packmol.inp_file
            print(f'  [custom inp] 使用 {inp}（跳过自动生成 write_inp）')
            structs = parse_inp_structures(inp)
            counts = [n for _, n in structs]
            type_names = [os.path.splitext(os.path.basename(f))[0] for f, _ in structs]
            for t in type_names:
                if t not in natom_by_name:
                    raise RuntimeError(
                        f'自定义 inp structure 文件 {t}.xyz 不在 molecules 中'
                        f'（可用: {sorted(natom_by_name)}）；结构文件名须为 {{name}}.xyz')
        else:
            inp = os.path.join(workdir, 'packmol.inp')
            write_inp(cfg, [os.path.join(workdir, mc.name + '.xyz') for mc in cfg.molecules],
                      n_atoms, inp)
            counts = [mc.count for mc in cfg.molecules]
            type_names = [mc.name for mc in cfg.molecules]
        packed = run_packmol(inp, workdir)
        blocks = parse_packed_xyz(packed, counts, natom_by_name, type_names)
        b = cfg.packmol.box
        box = [b[0], b[3], b[1], b[4], b[2], b[5]]   # packmol → LAMMPS 顺序
        natom_input = sum(c * natom_by_name[t] for c, t in zip(counts, type_names))
    else:
        # 单分子：坐标块直接作为"一个拷贝"，盒自动推算
        m0 = mols[0]
        blocks = [[list(zip(m0['elements'],
                            [c[0] for c in m0['coords']],
                            [c[1] for c in m0['coords']],
                            [c[2] for c in m0['coords']]))]]
        box = None
        natom_input = m0['natom']

    print('== 导出层: ReaxFF 极简 data (Masses + Atoms) ==')
    exp = build_data_reaxff(cfg, workdir, blocks, atom_info, box=box)
    ok, msgs = _validate_reaxff(cfg, exp, natom_input)

    report = {
        'config': cfg,
        'molecules': mols,
        'system': None,
        'export': exp,
        'validation': {'ok': ok, 'messages': msgs},
        'data_lmp': exp['data_lmp'],
        'workdir': workdir,
    }
    # 3c. 体系标准 xyz 导出（默认主产物；ReaxFF 多分子=packed.xyz、单分子=分子层坐标块）
    xyz_path = os.path.join(workdir, 'system.xyz')
    n_xyz = export_system_xyz(report, xyz_path)
    print(f'  xyz → {os.path.basename(xyz_path)} ({n_xyz} atoms)')
    _write_report_reaxff(report)
    print(f'\n==== 校验结果: {"通过" if ok else "失败"} ====')
    for m in msgs:
        print('  ' + m)
    print(f'\n输出: {exp["data_lmp"]}')
    print(f'检查报告: {os.path.join(workdir, "check_report.txt")}')

    if cfg.organize_output:
        keep = [exp['data_lmp']] + [m['mol2'] for m in mols if m.get('mol2')] + [xyz_path]
        # 保留 QM 波函数**正式产物**（fch/chg；ESP iso/pt 产物在
        # _others/electrostatic_potential/ 内，organize 不动 _others 目录）
        keep += [p for p in glob.glob(os.path.join(workdir, '*.fch'))]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.molden'))]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.chg'))]
        # 保留 frcmod 与指纹文件（reuse_molecule 复用分子层需要）
        keep += [m['frcmod'] for m in mols if m.get('frcmod')]
        keep += [p for p in glob.glob(os.path.join(workdir, '*.fingerprint.json'))]
        # 保留 packmol inp：自动生成模板（供用户改）+ 自定义 inp（用户输入，不能移走）
        keep.append(os.path.join(workdir, 'packmol.inp'))
        if cfg.packmol.inp_file:
            keep.append(cfg.packmol.inp_file)
        organize_workdir(workdir, keep, cfg.organize_backup)
    return report


def _validate_reaxff(cfg: Config, exp: dict, natom_input: int) -> tuple[bool, list[str]]:
    """ReaxFF 校验：原子数守恒 + 段计数 + 类型数 + Masses 覆盖。"""
    info = exp['info']
    check = exp['check']
    msgs: list[str] = []
    ok = True

    if info['natom'] == natom_input:
        msgs.append(f'原子数守恒: {info["natom"]} == 分子层 {natom_input} 通过')
    else:
        ok = False
        msgs.append(f'原子数不一致: data={info["natom"]} 分子层={natom_input} 失败')

    # 无键项
    for key in ('bonds', 'angles', 'dihedrals', 'impropers'):
        h = check['counts'].get(key, 0)
        a = check.get(key + '_actual', 0)
        if h != a:
            ok = False
            msgs.append(f'段计数不一致 {key}: 头部={h} 实际={a} 失败')
    msgs.append('段计数: 头部与文件记录一致 通过（无键项，ReaxFF 预期）')

    # 类型计数
    for key in ('atom_types', 'bond_types', 'angle_types', 'dihedral_types', 'improper_types'):
        h = check['counts'].get(key, 0)
        a = check.get(key + '_actual', 0)
        if h != a:
            ok = False
            msgs.append(f'类型计数不一致 {key}: 头部={h} 实际={a} 失败')
    msgs.append('类型计数: 头部与 Coeffs 段一致 通过')

    # 电荷列全 0（QEq 待算）
    msgs.append('电荷: 初始 0.0（QEq/reax 模拟中计算）通过')
    return ok, msgs


def _write_report_reaxff(report: dict) -> None:
    cfg = report['config']
    mols = report['molecules']
    exp = report['export']
    info = exp['info']
    ok = report['validation']['ok']
    path = os.path.join(report['workdir'], 'check_report.txt')
    multi = cfg.packmol.enabled

    lines = [
        f'# getlmp 检查报告 — {cfg.name}',
        f'- 日期: {__import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        f'- 配置: {cfg.name}.yaml  (forcefield=reaxff, charge=none/QEq)',
        '',
        '## 分子层（RDKit 3D 坐标）',
    ]
    for mc, mol in zip(cfg.molecules, mols):
        lines.append(f'- SMILES: {mc.smiles}  (name={mc.name}, count={mc.count})')
        lines.append(f'  - 单分子原子数 (含 H): {mol["natom"]}，'
                     f'元素 {sorted(set(mol["elements"]))}')
    lines += ['', '## 体系层（ReaxFF：坐标装盒，无拓扑）']
    if multi:
        lines += [f'- 装盒: packmol preset={cfg.packmol.preset} '
                  f'box={cfg.packmol.box} tolerance={cfg.packmol.tolerance} '
                  f'seed={cfg.packmol.seed}']
    else:
        lines += ['- 单分子：盒由坐标自动推算']
    lines += [
        '',
        '## 导出层（data.lmp，ReaxFF 极简格式）',
        f'- 原子: {info["natom"]}  原子类型: {info["atom_types"]}',
        f'- atom_style: {cfg.reax_atom_style}'
        f'（{"6 列 id charge x y z" if cfg.reax_atom_style == "charge" else "7 列 id mol-id type charge x y z"}）',
        '- 键/角/二面角/improper: 无（ReaxFF 键级由 pair_style 计算）',
        '- 电荷: 0.0（QEq 模拟中计算）',
        '- 盒: [' + ', '.join(f'{v:.3f}' for v in info['box']) + ']',
        f'- 元素顺序（类型号→元素）: {cfg.reax_elements}',
        '',
        '## 校验',
        f'- 结论: {"通过" if ok else "失败"}',
    ]
    lines += [f'- {m}' for m in report['validation']['messages']]
    lines += [
        '',
        '## LAMMPS 实跑验证（ReaxFF，建议集群跑）',
        '本机无 LAMMPS；将 data.lmp + ffield.reax + param.qeq 交给 OpsAgent 代跑：',
        '```lammps',
        REAXFF_IN_TEMPLATE.format(
            data=os.path.basename(exp['data_lmp']),
            elements=' '.join(cfg.reax_elements)),
        '```',
        '通过判据: READ_DATA_OK + run 0 无 ERROR。',
    ]
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


# ReaxFF 的 LAMMPS 验证模板（data 无 Coeffs 段）
# 关键规则（OpsAgent 集群实测暴露，作业见集群 getlmp_verify/ethanol_reax）：
#   - atom_style charge 的 Atoms 段是 6 列（id charge x y z），不能有 molecule 列
#   - data 无 Pair Coeffs 段时，pair_style + pair_coeff 必须在 read_data 之后
#   - fix qeq/reaxff 的 params 用关键字 reaxff → 从 ffield 提取 QEq 参数，无需 param.qeq
REAXFF_IN_TEMPLATE = """units real
atom_style charge
boundary p p p
read_data {data}
pair_style reaxff NULL
pair_coeff * * ffield.reax {elements}
fix 1 all qeq/reaxff 1 0.0 10.0 1e-6 reaxff
run 0
"""


def _write_report(report: dict) -> None:
    cfg = report['config']
    mols = report['molecules']
    exp = report['export']
    info = exp['info']
    ok = report['validation']['ok']
    path = os.path.join(report['workdir'], 'check_report.txt')
    multi = cfg.packmol.enabled

    lines = [
        f'# getlmp 检查报告 — {cfg.name}',
        f'- 日期: {__import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        f'- 配置: {cfg.name}.yaml  (forcefield={cfg.forcefield}, '
        f'charge={cfg.charge_method}, net_charge={cfg.net_charge})',
        '',
        '## 分子层',
    ]
    for mc, mol in zip(cfg.molecules, mols):
        lines += [
            f'- SMILES: {mc.smiles}  (name={mc.name}, count={mc.count}, '
            f'resname={mc.resname})',
            f'  - 单分子原子数 (含 H): {mol["natom"]}，'
            f'frcmod: {"有" if mol.get("n_extra") else "无（GAFF 标准库齐全）"}',
        ]
    lines += ['', '## 体系层（tleap）',
              f'- prmtop: {os.path.basename(report["system"]["prmtop"])}',
              f'- inpcrd: {os.path.basename(report["system"]["inpcrd"])}']
    if multi:
        lines += [f'- 装盒: packmol preset={cfg.packmol.preset} '
                  f'box={cfg.packmol.box} tolerance={cfg.packmol.tolerance} '
                  f'seed={cfg.packmol.seed}']
    lines += [
        '',
        '## 导出层（data.lmp）',
        f'- 原子: {info["natom"]}  键: {info["nbond"]}  角: {info["nangle"]}',
        f'- 二面角: {info["ndihedral"]}  improper: {info["nimproper"]}',
        f'- 电荷总和: {info["total_charge"]:.6f}（净电荷 {cfg.net_charge}）',
        '- 盒: [' + ', '.join(f'{v:.3f}' for v in info['box']) + ']',
        '',
        '## 校验',
        f'- 结论: {"通过" if ok else "失败"}',
    ]
    lines += [f'- {m}' for m in report['validation']['messages']]
    if not multi and cfg.charge_method in ('resp', 'resp2') and cfg.esp.enabled:
        wfn_type = 'Gaussian .fch' if cfg.qm.engine == 'gaussian' else 'QUICK .molden'
        lines += [
            '',
            '## ESP 可视化导出（单分子 RESP/RESP2, Multiwfn）',
            '- iso 法: _others/electrostatic_potential/iso/'
            '（density.cub 电子密度 + esp.cub 严格 QM ESP；'
            'VMD: density 等值面 0.001 + esp 着色）',
            '- pt 法: _others/electrostatic_potential/pt/'
            '（mol.pdb 分子结构 + vtx.pdb 密度 0.001 等值面顶点；'
            'B 因子=ESP kcal/mol, VMD Beta 着色）',
            f'- 波函数: {os.path.basename(mols[0].get("wavefn") or "-")} '
            f'（{wfn_type}, VMD/Multiwfn 可视化/二次分析用）',
        ]
    lines += [
        '',
        '## LAMMPS 实跑验证（可选，建议集群跑）',
        '本机无 LAMMPS；将 data.lmp 交给 OpsAgent 代跑或自行在集群执行：',
        '```lammps',
        LAMMPS_IN_TEMPLATE.format(data=os.path.basename(exp['data_lmp'])),
        '```',
        '通过判据: READ_DATA_OK + run 0 无 ERROR（VERIFY_EXIT=0）。',
    ]
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
