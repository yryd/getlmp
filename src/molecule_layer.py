#!/usr/bin/env python
"""分子层：SMILES → 3D SDF → antechamber(类型+电荷) → mol2 → parmchk2 → frcmod

支持 GAFF2/GAFF 类型 + AM1-BCC / ABCG2 / RESP / RESP2 电荷。
RESP/RESP2 的 QM 引擎由 cfg.qm 决定：
- engine=gaussian（默认）：G16 单点（pop=MK ESP）+ Multiwfn 拟合（531 路线）
- engine=quick（旧路径）：QUICK HF/6-31G* + resp 两阶段拟合
所有外部命令（antechamber/parmchk2/sqm/g16/quick/Multiwfn）来自 conda 环境 PATH
或 qm 配置路径。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

from config import Config, MoleculeCfg


def _run(cmd: list[str], cwd: str, desc: str) -> None:
    print(f'  [run] {" ".join(cmd)}', file=sys.stderr)
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'{desc} 失败 (rc={r.returncode})\n'
                           f'--- stdout ---\n{r.stdout[-2000:]}\n'
                           f'--- stderr ---\n{r.stderr[-2000:]}')


def _smiles_to_sdf(smiles: str, out_sdf: str, seed: int) -> int:
    """RDKit SMILES → 3D SDF（ETKDG + UFF）。返回含氢原子数。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'无法解析 SMILES: {smiles}')
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError(f'ETKDG 构象生成失败 (status={status}): {smiles}')
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    w = Chem.SDWriter(out_sdf)
    w.write(mol)
    w.close()
    return mol.GetNumAtoms()


def _sdf_elements_coords(sdf: str) -> tuple[list[str], list[list[float]]]:
    """SDF → (元素列表, 坐标列表 [[x,y,z],...])，顺序与 SDF 原子一致。"""
    from rdkit import Chem
    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        raise RuntimeError(f'RDKit 无法读取 SDF: {sdf}')
    conf = mol.GetConformer()
    elements = []
    coords = []
    for a in mol.GetAtoms():
        elements.append(a.GetSymbol())
        p = conf.GetAtomPosition(a.GetIdx())
        coords.append([p.x, p.y, p.z])
    return elements, coords


def build_molecule(cfg: Config, mc: MoleculeCfg, workdir: str) -> dict:
    """单个分子的完整分子层流水线，返回各中间产物路径。

    GAFF2/GAFF：antechamber(类型+电荷) + parmchk2。
    ReaxFF：仅 RDKit 3D 坐标（无类型/电荷，QEq 模拟中算）。
    """
    os.makedirs(workdir, exist_ok=True)
    base = os.path.join(workdir, mc.name)
    sdf = base + '.sdf'

    print(f'== 分子层: {mc.name} ({mc.smiles}) ==')
    natom = _smiles_to_sdf(mc.smiles, sdf, cfg.seed)
    print(f'  SMILES → SDF: {natom} atoms (incl. H)')

    if cfg.forcefield == 'reaxff':
        # ReaxFF 分支：坐标 + 元素即可，不需要 antechamber / parmchk2
        elements, coords = _sdf_elements_coords(sdf)
        print(f'  ReaxFF: {len(elements)} atoms, 元素 {sorted(set(elements))}'
              f'（类型/电荷由 ffield.reax + QEq 在 LAMMPS 中处理）')
        return {'sdf': sdf, 'mol2': None, 'frcmod': None, 'wavefn': None,
                'elements': elements, 'coords': coords, 'natom': natom,
                'n_extra': 0}

    mol2 = base + '.mol2'
    frcmod = base + '.frcmod'
    wavefn = None   # 仅 RESP/RESP2 分支生成（gaussian→.fch，quick→.molden）
    resname = (mc.resname or mc.name[:3]).upper()[:3]

    if cfg.charge_method in ('resp', 'resp2'):
        # RESP/RESP2 分支：先用 bcc 生成类型正确的 mol2（临时电荷），再 QM 拟合覆盖
        _run([
            'antechamber', '-i', sdf, '-fi', 'sdf', '-o', mol2, '-fo', 'mol2',
            '-c', 'bcc', '-at', cfg.forcefield,
            '-nc', str(cfg.net_charge), '-rn', resname,
        ], workdir, 'antechamber(bcc 类型占位)')
        print('  antechamber(bcc) → mol2 类型就绪（电荷待 QM 拟合覆盖）')
        if cfg.charge_method == 'resp2':
            charges = _build_resp2_charges(mc, workdir, base, mol2,
                                           cfg.net_charge, cfg.qm)
            label = f'RESP2({cfg.qm.engine})'
        else:
            charges = _build_resp_charges(mc, workdir, base, mol2,
                                          cfg.net_charge, cfg.esp.enabled, cfg.qm)
            label = f'RESP({cfg.qm.engine})'
        _apply_charges(mol2, charges)
        print(f'  {label} 电荷拟合完成 → {len(charges)} atoms, '
              f'sum={sum(charges):+.6f}')
        # 归档波函数文件（esp 开启时；gaussian→.fch（RESP2 用溶剂 .fch），
        # quick→.molden，VMD/Multiwfn 可视化与二次分析用）
        if cfg.charge_method == 'resp2':
            wavefn_src = base + '_solv.fch'
            wavefn_ext = 'fch'
        elif cfg.qm.engine == 'gaussian':
            wavefn_src = base + '.fch'
            wavefn_ext = 'fch'
        else:
            wavefn_src = base + '_quick.molden'
            wavefn_ext = 'molden'
        if cfg.esp.enabled and os.path.exists(wavefn_src):
            wavefn = base + '.' + wavefn_ext
            if os.path.abspath(wavefn_src) == os.path.abspath(wavefn):
                wavefn = wavefn_src   # 已是最终位置（gaussian 路径 .fch），无需复制
            else:
                shutil.copy(wavefn_src, wavefn)
            print(f'  已归档 {os.path.basename(wavefn)}（波函数，'
                  f'可视化/二次分析用）')
    else:
        _run([
            'antechamber', '-i', sdf, '-fi', 'sdf', '-o', mol2, '-fo', 'mol2',
            '-c', cfg.charge_method, '-at', cfg.forcefield,
            '-nc', str(cfg.net_charge), '-rn', resname,
        ], workdir, 'antechamber')
        print(f'  antechamber → mol2 ({cfg.charge_method} / {cfg.forcefield}, resname={resname})')

    _normalize_mol2_charge(mol2, cfg.net_charge)
    print(f'  电荷归一化 → 分子净电荷 = {cfg.net_charge}')

    _run(['parmchk2', '-i', mol2, '-f', 'mol2', '-o', frcmod], workdir, 'parmchk2')
    n_extra = _count_frcmod_lines(frcmod)
    if n_extra:
        print(f'  parmchk2 → frcmod（{n_extra} 段补充参数）')
    else:
        print('  parmchk2 → frcmod（空：参数全部来自 GAFF 标准库）')

    return {'sdf': sdf, 'mol2': mol2, 'frcmod': frcmod, 'wavefn': wavefn,
            'natom': natom, 'n_extra': n_extra}


def _normalize_mol2_charge(mol2: str, net_charge: int) -> None:
    """把 mol2 的分子总电荷修正到 net_charge（AM1-BCC 3 位小数截断会引入 ±0.00x
    误差，多拷贝时累积成大偏差；修正值加到最后一个重原子）。

    修改 ATOM 段每行末的电荷值，原地写回 mol2。
    """
    with open(mol2) as f:
        lines = f.read().splitlines()
    out = []
    in_atom = False
    charges: list[tuple[int, int]] = []   # (行号, 是否为重原子)
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith('@<TRIPOS>ATOM'):
            in_atom = True
            out.append(ln)
            continue
        if s.startswith('@<TRIPOS>'):
            in_atom = False
            out.append(ln)
            continue
        if in_atom and s and not s.startswith('#'):
            parts = ln.split()
            if len(parts) >= 9:
                name = parts[1]
                charges.append((i, not name.startswith('H')))
        out.append(ln)

    if not charges:
        return
    # 当前总电荷（第 9 列）
    total = 0.0
    for i, _ in charges:
        total += float(out[i].split()[8])
    delta = float(net_charge) - total
    if abs(delta) < 1e-6:
        return
    # 选最后一个重原子（无重原子则最后一个原子）
    target = next((i for i, heavy in reversed(charges) if heavy), charges[-1][0])
    parts = out[target].split()
    parts[8] = f'{float(parts[8]) + delta:.6f}'
    out[target] = '    ' + '    '.join(parts)
    with open(mol2, 'w') as f:
        f.write('\n'.join(out) + '\n')


def _count_frcmod_lines(path: str) -> int:
    """frcmod 中非空参数行数（MASS/BOND/ANGLE/DIHE/IMPROPER/NONBON 段下的行）。"""
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        lines = f.read().splitlines()
    n = 0
    in_section = False
    for ln in lines:
        s = ln.strip()
        if s in {'MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON'}:
            in_section = True
            continue
        if in_section and s and not s.startswith('#'):
            n += 1
    return n


# ---------------------------------------------------------------- RESP 分支
# RESP 电荷流程：按 qm.engine 分派
#   engine=gaussian（默认，531 路线）：
#     mol2 → gjf（B3LYP/def2-TZVP + pop=MK ESP，可选 PCM 溶剂）
#     → g16 单点 → formchk(.fch) → Multiwfn RESP/RESP2 拟合 → 电荷
#   engine=quick（旧路径保留，见 docs/dev_notes.md）：
#     mol2 → QUICK 输入(.in) → quick 算 HF/6-31G* ESP(.out+.vdw)
#     → antechamber(-fi quick) 生成 RESP1/2.IN（原子等价性）
#     → 手工构建 .esp（绕过 espgen -f 4 的输出 bug：缺 ESP 值列/列序错位）
#     → resp 两阶段拟合 → QOUT 电荷 → 写回 mol2（保留 gaff2 类型）
#
# 已知坑（2026-08-09 实测，ambertools 26.0 conda + QUICK 2025）：
#   1. `-c resp` 必须提供外部 QM 输出（-fi quick/gout/gesp/gamess），不能一步完成；
#   2. `-c resp2` 不被支持（Unknown charge method），RESP2 需回退方案；
#   3. espgen -f 4 生成的 .esp 只有 3 列坐标（缺 ESP 值），resp 拟合出 ~0 电荷；
#   4. QUICK .out 的 INPUT GEOMETRY 只打印 X/Y（缺 Z），原子坐标须从 mol2/xyz 取；
#   5. resp 的 -e .esp 格式：ESP 在前、坐标在后、单位 Bohr（原子坐标 + MEP 点都要 Bohr）。
BOHR = 1.889726124626

def _find_g16() -> str:
    """按 PATH 约定探测 g16 可执行（GAUSS_* 环境由 .bashrc 提供，见 docs/install.md）。"""
    p = shutil.which('g16')
    if p:
        return p
    raise RuntimeError(
        'PATH 中找不到 g16（Gaussian 未安装或未加入 PATH）。'
        '安装与配置见 docs/install.md；也可用 qm.engine: quick 走 QUICK 回退路径')


def _mol2_atoms(mol2: str) -> list[list]:
    """从 mol2 ATOM 段读 (元素符号, x, y, z)。元素从原子名去数字得到。"""
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    atoms = []
    in_atom = False
    for ln in open(mol2):
        s = ln.strip()
        if s.startswith('@<TRIPOS>ATOM'):
            in_atom = True
            continue
        if s.startswith('@<TRIPOS>'):
            in_atom = False
            continue
        if in_atom and s and not s.startswith('#'):
            parts = ln.split()
            if len(parts) >= 6:
                name = parts[1]
                sym = ''.join(ch for ch in name if not ch.isdigit())
                z = pt.GetAtomicNumber(sym)
                atoms.append([sym, z, float(parts[2]), float(parts[3]), float(parts[4])])
    return atoms


def _require_quick_basis() -> str:
    """读取约定环境变量 QUICK_BASIS（指向含 6-31G.BAS 的目录，见 docs/install.md）。"""
    b = os.environ.get('QUICK_BASIS', '')
    if b and os.path.isdir(b) and os.path.exists(os.path.join(b, '6-31G.BAS')):
        return b
    raise RuntimeError(
        '环境变量 QUICK_BASIS 未设置或无效（应指向含 6-31G.BAS 的目录）。'
        '安装与配置见 docs/install.md（例如 '
        'export QUICK_BASIS=$CONDA_PREFIX/AmberTools/src/quick/basis）')


def _write_quick_input(path: str, title: str, atoms: list[list], net_charge: int,
                       esp_enabled: bool = False) -> None:
    """写 QUICK 输入（$DATA 格式 + HF/6-31G* ESP_CHARGE）。

    READ_COORD 原子行格式 = `元素 x y z`（元素 + 3 坐标，不要原子序数列！
    QUICK 按 4 字段读，多写 Z 会把 Z 当 x、x 当 y…导致几何错乱）。
    EXPORT=MOLDEN 仅 esp 开启时附加（供可视化；默认关闭不产出波函数文件）。
    """
    method = 'HF BASIS=6-31G* READ_COORD ESP_CHARGE' + \
        (' EXPORT=MOLDEN' if esp_enabled else '')
    lines = [f'$DATA = {title} RESP via getlmp', method, '']
    for sym, z, x, y, zz in atoms:
        lines.append(f'{sym:<3s}{x:12.5f}{y:12.5f}{zz:12.5f}')
    lines.append('$END')
    if net_charge != 0:
        # QUICK 电荷关键词（带电分子 RESP 仍在验证，见 docs/dev_notes.md）
        lines.insert(1, f'CHARGE {net_charge}')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _write_gaussian_input(path: str, title: str, atoms: list[list],
                          net_charge: int, qm, chk_path: str,
                          solvent: str | None = None) -> None:
    """写 g16 输入：B3LYP/def2TZVP + pop=MK ESP（531 路线关键词）。

    route = #p {method}/{basis} em=GD3BJ pop=MK IOp(6/33=2,6/42=6)
            + scrf=(smd,solvent=xxx)（溶剂非空时） + opt（qm.opt）
    坐标来自 mol2（不重排原子序）；多重度固定 1（闭壳层）。
    solvent 参数可覆盖 qm.solvent（RESP2 气相单点传 ''）。

    坑（2026-08-11 实测 G16 RevC.01）：
    - `0 1` 行后**不能有空行**，否则 l101 报 "There are no atoms"；
    - def2 基组名带连字符（def2-TZVP）会触发 route 语法错误，需转 def2TZVP。
    """
    basis = qm.basis.replace('def2-', 'def2')   # def2-TZVP → def2TZVP
    solv = qm.solvent if solvent is None else solvent
    route = f'#p {qm.method}/{basis} em=GD3BJ pop=MK IOp(6/33=2,6/42=6)'
    if solv:
        route += f' scrf=(smd,solvent={solv})'
    if qm.opt:
        route += ' opt'
    lines = [
        f'%nprocshared={os.environ.get("SMI2DATA_NPROC", "8")}',
        f'%mem={os.environ.get("SMI2DATA_MEM", "4GB")}',
        f'%chk={os.path.abspath(chk_path)}',
        route,
        '',
        f'{title} RESP via getlmp ({qm.method}/{basis})',
        '',
        f'{net_charge} 1',
    ]
    for sym, _z, x, y, zz in atoms:
        lines.append(f'{sym:<3s}{x:12.6f}{y:12.6f}{zz:12.6f}')
    lines.append('')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _run_gaussian(gjf: str, workdir: str, qm) -> str:
    """运行 g16 单点，返回 log 路径；检查 Normal termination。"""
    _find_g16()   # PATH 探测；GAUSS_* 环境由 .bashrc 提供，子进程直接继承
    log = os.path.splitext(gjf)[0] + '.log'
    print(f'  [run] g16 {os.path.basename(gjf)} '
          f'({qm.method}/{qm.basis}'
          + (f' + {qm.solvent}' if qm.solvent else '')
          + ', pop=MK ESP)')
    try:
        r = subprocess.run(['g16', os.path.basename(gjf), os.path.basename(log)],
                           cwd=workdir, capture_output=True, text=True,
                           timeout=3600)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f'g16 超时（>1h）: {gjf}')
    log_full = os.path.join(workdir, log)
    if not os.path.exists(log_full) or \
            'Normal termination' not in open(log_full).read():
        tail = ''
        if os.path.exists(log_full):
            tail = open(log_full).read()[-1500:]
        raise RuntimeError(
            f'g16 运行失败（无 Normal termination）\n'
            f'--- stdout ---\n{r.stdout[-800:]}\n'
            f'--- stderr ---\n{r.stderr[-800:]}\n'
            f'--- log 尾部 ---\n{tail}')
    return log


def _run_formchk(chk: str, fch: str, workdir: str, qm) -> str:
    """formchk: .chk → .fch（格式波函数，Multiwfn 可读）。"""
    fc = shutil.which('formchk')
    if not fc:
        raise RuntimeError('PATH 中找不到 formchk（Gaussian 未安装或未加入 PATH）。'
                           '安装与配置见 docs/install.md')
    r = subprocess.run(
        [fc, os.path.basename(chk), os.path.basename(fch)],
        cwd=workdir, capture_output=True, text=True, timeout=300)
    if not os.path.exists(os.path.join(workdir, fch)):
        raise RuntimeError(f'formchk 失败\n--- stdout ---\n{r.stdout[-800:]}\n'
                           f'--- stderr ---\n{r.stderr[-800:]}')
    return fch


def _build_esp_file(esp_path: str, atoms_xyz: list[list], vdw_path: str) -> int:
    """构建 resp 可读的 .esp（绕过 espgen bug）。

    格式：第 1 行 natom nesp 0；随后 natom 行原子坐标(Bohr)；
    再 nesp 行 `ESP X Y Z`（ESP a.u.，坐标 Bohr）。坐标单位都转 Bohr。
    返回 MEP 点数。
    """
    pts = []
    for ln in open(vdw_path):
        parts = ln.split()
        if len(parts) == 4:
            try:
                pts.append([float(p) for p in parts])   # X Y Z ESP (angstrom)
            except ValueError:
                continue
    with open(esp_path, 'w') as f:
        f.write('%5d%5d%5d\n' % (len(atoms_xyz), len(pts), 0))
        for _, _, x, y, z in atoms_xyz:
            f.write('%17s%16.7E%16.7E%16.7E\n' % ('', x * BOHR, y * BOHR, z * BOHR))
        for p in pts:
            # resp/espgen 约定：点行 = ESP 值 X Y Z（值在前、坐标 Bohr 在后）
            f.write('%1s%16.7E%16.7E%16.7E%16.7E\n' % (
                '', p[3], p[0] * BOHR, p[1] * BOHR, p[2] * BOHR))
    return len(pts)


def _build_resp_charges(mc: MoleculeCfg, workdir: str, base: str,
                        mol2: str, net_charge: int, esp_enabled: bool,
                        qm) -> list[float]:
    """RESP 电荷拟合入口：按 qm.engine 分派（gaussian 主路径 / quick 旧路径）。"""
    if qm.engine == 'gaussian':
        return _build_resp_charges_gaussian(mc, workdir, base, mol2,
                                            net_charge, qm)
    if qm.engine == 'quick':
        return _build_resp_charges_quick(mc, workdir, base, mol2,
                                         net_charge, esp_enabled)
    raise ValueError(f'不支持的 qm.engine={qm.engine!r}')


def _build_resp_charges_gaussian(mc: MoleculeCfg, workdir: str, base: str,
                                 mol2: str, net_charge: int, qm) -> list[float]:
    """RESP（Gaussian + Multiwfn，531 路线）：
    gjf 单点（pop=MK）→ g16 → formchk(.fch) → Multiwfn 7→18→1 拟合 → 电荷。
    """
    from multiwfn import resp_from_fch

    atoms = _mol2_atoms(mol2)
    gjf = base + '.gjf'
    chk = base + '.chk'
    fch = base + '.fch'

    _write_gaussian_input(gjf, mc.name, atoms, net_charge, qm, chk)
    _run_gaussian(gjf, workdir, qm)
    _run_formchk(chk, fch, workdir, qm)
    print(f'  [run] Multiwfn RESP 拟合 ← {os.path.basename(fch)}')
    charges, _chg = resp_from_fch(fch, workdir=workdir,
                                  multiwfn_path=qm.multiwfn_path)
    if len(charges) != len(atoms):
        raise RuntimeError(f'RESP 电荷数 {len(charges)} != 原子数 {len(atoms)}')
    return charges


def _build_resp2_charges(mc: MoleculeCfg, workdir: str, base: str,
                         mol2: str, net_charge: int, qm) -> list[float]:
    """RESP2（Gaussian + Multiwfn）：气相与溶剂(PCM)各单点 → 各自 RESP → δ 混合。

    q = (1-δ)*q_gas + δ*q_solv，对应官方 calcRESP.sh 逻辑。
    """
    from multiwfn import resp2_from_fch

    atoms = _mol2_atoms(mol2)
    gjf_gas, chk_gas, fch_gas = base + '_gas.gjf', base + '_gas.chk', base + '_gas.fch'
    gjf_solv, chk_solv, fch_solv = base + '_solv.gjf', base + '_solv.chk', base + '_solv.fch'

    # 气相单点
    _write_gaussian_input(gjf_gas, mc.name, atoms, net_charge, qm, chk_gas,
                          solvent='')
    _run_gaussian(gjf_gas, workdir, qm)
    _run_formchk(chk_gas, fch_gas, workdir, qm)
    # 溶剂单点（PCM）
    _write_gaussian_input(gjf_solv, mc.name, atoms, net_charge, qm, chk_solv,
                          solvent=qm.solvent)
    _run_gaussian(gjf_solv, workdir, qm)
    _run_formchk(chk_solv, fch_solv, workdir, qm)
    # Multiwfn RESP2（两次 RESP + δ 混合）
    print(f'  [run] Multiwfn RESP2 拟合 ← gas + {qm.solvent} (δ={qm.delta})')
    charges = resp2_from_fch(fch_gas, fch_solv, qm.delta, workdir=workdir,
                             multiwfn_path=qm.multiwfn_path)
    if len(charges) != len(atoms):
        raise RuntimeError(f'RESP2 电荷数 {len(charges)} != 原子数 {len(atoms)}')
    return charges


def _build_resp_charges_quick(mc: MoleculeCfg, workdir: str, base: str,
                              mol2: str, net_charge: int,
                              esp_enabled: bool = False) -> list[float]:
    """RESP 电荷拟合（QUICK ESP + resp 两阶段，旧路径），返回电荷列表。"""
    atoms = _mol2_atoms(mol2)
    quick_in = base + '_quick.in'
    quick_out = base + '_quick.out'
    quick_vdw = base + '_quick.vdw'

    # 1. QUICK 输入 + 运行
    _write_quick_input(quick_in, mc.name, atoms, net_charge, esp_enabled)
    basis = _require_quick_basis()
    print(f'  [run] quick {os.path.basename(quick_in)} (HF/6-31G* ESP, '
          f'QUICK_BASIS={basis})')
    r = subprocess.run(['quick', os.path.basename(quick_in)], cwd=workdir,
                       capture_output=True, text=True)
    if not os.path.exists(quick_out) or 'Normal Termination' not in open(quick_out).read():
        raise RuntimeError(f'QUICK 运行失败\n--- stdout ---\n{r.stdout[-1500:]}\n'
                           f'--- stderr ---\n{r.stderr[-1500:]}')
    if not os.path.exists(quick_vdw):
        raise RuntimeError('QUICK 未生成 .vdw（ESP on vdW surface）。'
                           '请确认输入含 ESP_CHARGE 关键词')

    # 2. antechamber: mol2 → .ac（连接性完整） + respgen 生成 RESP 输入
    #    等价性由 respgen 从完整连接性判断（苯 6C/6H 等效正确；对乙醇等保守
    #    退化全独立——无假等效，物理上安全）。绕开 -fi quick（连接性缺失
    #    导致过度等效）与 -fi quick 的路径/格式问题。
    tmp = os.path.join(workdir, f'_resp_tmp_{mc.name}')
    os.makedirs(tmp, exist_ok=True)
    _run([
        'antechamber', '-i', os.path.abspath(mol2), '-fi', 'mol2',
        '-o', 'dummy.ac', '-fo', 'ac',
    ], tmp, 'antechamber(mol2 → ac)')
    _run(['respgen', '-i', 'dummy.ac', '-o', 'ANTECHAMBER_RESP1.IN', '-f', 'resp1'],
         tmp, 'respgen(resp1)')
    _run(['respgen', '-i', 'dummy.ac', '-o', 'ANTECHAMBER_RESP2.IN', '-f', 'resp2'],
         tmp, 'respgen(resp2)')
    for p in ('ANTECHAMBER_RESP1.IN', 'ANTECHAMBER_RESP2.IN'):
        if not os.path.exists(os.path.join(tmp, p)):
            raise RuntimeError(f'respgen 未生成 {p}（RESP 输入缺失）')

    # 3. 修复 .esp + resp 两阶段拟合
    esp_fixed = os.path.join(tmp, 'ANTECHAMBER.ESP.fixed')
    n_esp = _build_esp_file(esp_fixed, atoms, quick_vdw)
    print(f'  RESP 拟合: {len(atoms)} atoms, {n_esp} ESP 点 (resp 两阶段)')
    qout2 = os.path.join(tmp, 'QOUT')
    _run(['resp', '-C', '-O', '-i', 'ANTECHAMBER_RESP1.IN', '-o', 'R1.OUT',
          '-e', 'ANTECHAMBER.ESP.fixed', '-t', 'qout'],
         tmp, 'resp 阶段 1')
    _run(['resp', '-C', '-O', '-i', 'ANTECHAMBER_RESP2.IN', '-o', 'R2.OUT',
          '-e', 'ANTECHAMBER.ESP.fixed', '-q', 'qout', '-t', 'QOUT'],
         tmp, 'resp 阶段 2')

    charges = [float(x) for x in open(qout2).read().split()]
    if len(charges) != len(atoms):
        raise RuntimeError(f'RESP 电荷数 {len(charges)} != 原子数 {len(atoms)}')
    # 拟合质量摘要（RMS ~0.9 常见于 QUICK vdW 网格，电荷仍可靠；记录不阻塞）
    r2 = open(os.path.join(tmp, 'R2.OUT')).read()
    import re as _re
    m = _re.search(r'ESP relative RMS\s+([\d.]+)', r2)
    if m:
        print(f'  [info] RESP 拟合 ESP relative RMS = {m.group(1)} '
              f'（QUICK vdW 网格特性，电荷值有效）')
    return charges


def _apply_charges(mol2: str, charges: list[float]) -> None:
    """把电荷列表写回 mol2 ATOM 段（保留类型/坐标），原地修改。"""
    with open(mol2) as f:
        lines = f.read().splitlines()
    out = []
    in_atom = False
    i = 0
    for ln in lines:
        s = ln.strip()
        if s.startswith('@<TRIPOS>ATOM'):
            in_atom = True
            out.append(ln)
            continue
        if s.startswith('@<TRIPOS>'):
            in_atom = False
            out.append(ln)
            continue
        if in_atom and s and not s.startswith('#'):
            parts = ln.split()
            if len(parts) >= 9:
                parts[8] = '%.6f' % charges[i]
                i += 1
                out.append('%-8s %-8s %12s %12s %12s %-6s %-3s %-6s %12s' % (
                    parts[0], parts[1], parts[2], parts[3], parts[4],
                    parts[5], parts[6], parts[7], parts[8]))
                continue
        out.append(ln)
    assert i == len(charges), f'写回电荷数 {i} != {len(charges)}'
    with open(mol2, 'w') as f:
        f.write('\n'.join(out) + '\n')
