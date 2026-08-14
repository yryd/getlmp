#!/usr/bin/env python
"""Multiwfn 封装层：RESP / RESP2 电荷拟合 + ESP 可视化产物（iso/pt 两套）。

交互序列（2026-08-11 实测 Multiwfn 2026.7.15 noGUI）：
- RESP 电荷：  主功能 7 → 18 → 1 → y → 0 → 0 → q      → {base}.chg
- RESP2：两次 RESP（气相+溶剂）后按 δ 混合（calcRESP.sh 逻辑）→ RESP2.chg
- pt 法：      主功能 12 → 0 → 5(mol.pdb) → 6(vtx.pdb) → mol.pdb + vtx.pdb
                （mol.pdb 为分子结构；vtx.pdb 为密度 0.001 闭合等值面顶点，
                  B 因子字段存 ESP，单位 kcal/mol，VMD Beta 着色）
- iso 法 ESP： 主功能 12 → 3(间距0.4Å) → 0 → 13 → mapfunc.cub（QM ESP 网格）
- iso 法密度： 主功能 12 → 3(间距0.4Å) → 2 → 11(ED) → 0 → 13 → mapfunc.cub
                （两次运行网格一致：同一表面格点空间，iso 渲染直接叠加）
                （2026-08-14 修复：Multiwfn 默认表面格点间距 0.25 Å，
                  对 ≥1000 基函数大分子表面点数爆炸 → 900s 超时；
                  插入 3→spacing=0.4 后点数约为 1/4，乙醇实测 2.3 秒）

波函数输入支持 .fch（Gaussian）与 .molden（QUICK），Multiwfn 按扩展名自动识别。

Multiwfn 交互式输入用 stdin 管道喂入；退出时可能报 Fortran I/O 错误
（exit 59，stdin 耗尽），只要目标产物文件已生成即视为成功。

已知坑（2026-08-12 实测）：
- 大分子（约 ≥1000 基函数，def2TZVP 含 F 轨道）fch 加载到
  "Generating density matrix based on SCF orbitals" 时，默认 8 MB 栈
  上限会段错误：forrtl: severe (174): SIGSEGV，退出码 174，不产出 .chg。
  _run_multiwfn 已通过 preexec_fn 解除子进程栈限制（ulimit -s unlimited），
  并设置 OMP_STACKSIZE/KMP_STACKSIZE=1G（手册 2.1.2）。若仍见 174，
  先确认 shell 里 `ulimit -s unlimited` 后再启动管道。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

try:
    import resource  # Unix 专用：解除子进程栈限制（Multiwfn 大分子 SIGSEGV 修复）
except ImportError:  # pragma: no cover - Windows 兜底
    resource = None

# Multiwfn RESP 标准输入（官方 calcRESP.sh）
_RESP_INPUT = '7\n18\n1\ny\n0\n0\nq\n'
# iso 法 ESP：12→3 设表面格点间距 → 0 开始分析（映射函数默认 ESP）→ 13 导出 mapfunc.cub
def _esp_cube_input(spacing: float = 0.4) -> str:
    """iso 法 ESP 交互序列：12→3 设间距→0 分析→13 导出 mapfunc.cub。

    spacing: 表面格点间距 Å（默认 0.4）。Multiwfn 默认 0.25，对 ≥1000 基函数
    大分子表面点数爆炸、QM ESP 积分极慢，900s 会超时；0.4 Å 点数约为 1/4
    （乙醇实测 36480 点 / 2.3 秒）。小分子可显式传 0.25 更精细。
    """
    return f'12\n3\n{spacing}\n0\n13\n-1\n0\nq\n'


# iso 法密度：12→3 设间距 → 2→11 改映射函数为 Electron density → 0 分析 → 13 导出
def _esp_density_cube_input(spacing: float = 0.4) -> str:
    """iso 法电子密度交互序列（与 esp_cube 同一表面格点空间，可直接叠加渲染）。"""
    return f'12\n3\n{spacing}\n2\n11\n0\n13\n-1\n0\nq\n'


# pt 法：12→3 设间距 → 0 分析 → 5 导出 mol.pdb（需显式文件名）→ 6 导出 vtx.pdb → -1 返回
def _esp_pt_input(spacing: float = 0.25) -> str:
    """pt 法交互序列：12→3 设间距→0 分析→5 mol.pdb→6 vtx.pdb→-1 返回。

    spacing 默认 0.25（pt 法默认关，用精细档；调用方通常传与 iso 相同的值）。
    """
    return f'12\n3\n{spacing}\n0\n5\nmol.pdb\n6\n-1\n0\nq\n'


def auto_esp_params(natom: int) -> tuple[float, int | None]:
    """按原子数自动选 iso 法 ESP 参数（2026-08-14 实测标定 + 外推）。

    | 原子数 | spacing | timeout |
    |--------|---------|---------|
    | ≤ 20   | 0.25 Å  | 600 s   |
    | 21-40  | 0.3 Å   | 1800 s  |
    | > 40   | 0.4 Å   | 3600 s  |

    依据：乙醇(9 原子) 0.25 实测 9s；MPC(60 原子) 0.25 超时(>900s)、
    0.4 约 18min、0.6 约 7min。中间档按
    耗时 ∝ 表面积/spacing³ × 基函数数² 外推（见 auto_esp_params 调用处注释）。
    timeout=None 表示不限。
    """
    if natom <= 20:
        return 0.25, 600
    if natom <= 40:
        return 0.3, 1800
    return 0.4, 3600


def find_multiwfn(cfg_path: str = '') -> str:
    """探测 Multiwfn_noGUI 可执行路径（PATH 约定，见 docs/install.md）。

    顺序：config 显式路径（qm.multiwfn_path）→ PATH。
    找不到时报错并给出提示。
    """
    if cfg_path:
        if os.path.isfile(cfg_path):
            return cfg_path
        raise FileNotFoundError(f'multiwfn_path 不存在: {cfg_path}')
    for name in ('Multiwfn_noGUI', 'multiwfn'):
        p = shutil.which(name)
        if p:
            return p
    raise FileNotFoundError(
        'PATH 中找不到 Multiwfn_noGUI（未安装或未加入 PATH）。'
        '安装与配置见 docs/install.md；也可在 yaml 用 qm.multiwfn_path 显式指定')


def _unlimited_stack() -> None:
    """preexec_fn：把子进程栈上限设为无限。

    Multiwfn 读大分子（约 ≥1000 基函数，含 F 轨道）fch 生成密度矩阵时，
    默认 8 MB 栈会段错误（forrtl severe 174 SIGSEGV）；ulimit -s unlimited
    后正常（Multiwfn 手册 2.1.2，论坛 id=207 官方答复）。
    """
    if resource is not None:
        resource.setrlimit(resource.RLIMIT_STACK,
                           (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


def _run_multiwfn(mwfn: str, wfn_file: str, stdin_text: str,
                  cwd: str, timeout: int | None = 3600) -> str:
    """以 stdin 管道喂输入运行 Multiwfn，返回 stdout 文本。

    退出码非 0（如 59 = stdin 耗尽）不视为失败——由调用方检查产物文件。
    timeout 默认 3600s：大分子 iso 网格 QM ESP 积分可能数十分钟（见
    _esp_cube_input 注释），900s 会误杀（2026-08-14 修复）；None=不限。
    """
    env = dict(os.environ)
    env.setdefault('OMP_STACKSIZE', '1G')   # 手册 2.1.2：OpenMP 线程栈
    env.setdefault('KMP_STACKSIZE', '1G')   # Intel 编译器版本同样适用
    r = subprocess.run(
        [mwfn, os.path.abspath(wfn_file), '-ispecial', '1'],
        cwd=cwd, input=stdin_text, capture_output=True, text=True,
        timeout=timeout, env=env, preexec_fn=_unlimited_stack,
    )
    if r.returncode != 0 and r.returncode != 59:
        print(f'  [warn] Multiwfn 退出码 {r.returncode}（可能仍已产出文件）',
              file=sys.stderr)
    return r.stdout


def _parse_chg(path: str) -> list[float]:
    """解析 Multiwfn .chg：每行 name X Y Z charge → [charge,...]。"""
    chgs = []
    with open(path) as f:
        for ln in f:
            parts = ln.split()
            if len(parts) >= 5:
                try:
                    chgs.append(float(parts[4]))
                except ValueError:
                    continue
    if not chgs:
        raise RuntimeError(f'Multiwfn .chg 未解析到电荷: {path}')
    return chgs


def resp_from_fch(fch: str, workdir: str = '.',
                  multiwfn_path: str = '') -> tuple[list[float], str]:
    """对含 MK ESP 数据的 .fch 拟合 RESP 电荷。

    返回 (charges, chg_path)。chg_path 为 Multiwfn 输出的 {base}.chg。
    """
    mwfn = find_multiwfn(multiwfn_path)
    base = os.path.splitext(os.path.basename(fch))[0]
    chg_path = os.path.join(workdir, base + '.chg')
    if os.path.exists(chg_path):
        os.remove(chg_path)
    _run_multiwfn(mwfn, fch, _RESP_INPUT, cwd=workdir)
    if not os.path.exists(chg_path):
        raise RuntimeError(f'Multiwfn RESP 未生成电荷文件: {chg_path}')
    return _parse_chg(chg_path), chg_path


def resp2_from_fch(fch_gas: str, fch_solv: str, delta: float = 0.5,
                   workdir: str = '.', multiwfn_path: str = '') -> list[float]:
    """RESP2：气相与溶剂各拟合 RESP，再按 q = (1-δ)*q_gas + δ*q_solv 混合。

    对应官方 calcRESP.sh 逻辑（Multiwfn 无单命令 RESP2，两次 RESP + awk）。
    """
    q_gas, _ = resp_from_fch(fch_gas, workdir, multiwfn_path)
    q_solv, _ = resp_from_fch(fch_solv, workdir, multiwfn_path)
    if len(q_gas) != len(q_solv):
        raise RuntimeError(f'RESP2 原子数不一致: gas={len(q_gas)} solv={len(q_solv)}')
    return [(1 - delta) * a + delta * b for a, b in zip(q_gas, q_solv)]


def esp_pt(wfn: str, mol_pdb_out: str, vtx_pdb_out: str,
           workdir: str = '.', multiwfn_path: str = '',
           spacing: float = 0.25, timeout: int | None = 3600) -> tuple[str, str]:
    """pt 法：一次 Multiwfn 运行导出分子结构 mol.pdb + ESP 曲面顶点 vtx.pdb。

    vtx.pdb：密度 0.001 闭合等值面顶点，B 因子 = ESP kcal/mol（VMD Beta 着色）。
    spacing: 表面格点间距 Å（默认 0.25；通常与 iso 法保持一致）。
    timeout: 单次 Multiwfn 调用超时秒；None=不限。
    返回 (mol_pdb_out, vtx_pdb_out)。
    """
    mwfn = find_multiwfn(multiwfn_path)
    for tmp in ('mol.pdb', 'vtx.pdb'):
        p = os.path.join(workdir, tmp)
        if os.path.exists(p):
            os.remove(p)
    _run_multiwfn(mwfn, wfn, _esp_pt_input(spacing), cwd=workdir, timeout=timeout)
    src_mol = os.path.join(workdir, 'mol.pdb')
    src_vtx = os.path.join(workdir, 'vtx.pdb')
    if not os.path.exists(src_mol) or not os.path.exists(src_vtx):
        raise RuntimeError(f'Multiwfn pt 法未同时生成 mol.pdb/vtx.pdb: '
                           f'{(src_mol, src_vtx)}')
    shutil.copy(src_mol, mol_pdb_out)
    shutil.copy(src_vtx, vtx_pdb_out)
    return mol_pdb_out, vtx_pdb_out


def esp_cube(wfn: str, out_cub: str, workdir: str = '.',
             multiwfn_path: str = '', spacing: float = 0.4,
             timeout: int | None = 3600) -> str:
    """iso 法 ESP 网格：严格 QM ESP cube（mapfunc.cub，QM 电子云积分）。

    spacing: 表面格点间距 Å（默认 0.4，大分子加速；小分子可设 0.25 更精细）。
    timeout: 单次 Multiwfn 调用超时秒；None=不限。
    返回输出路径。
    """
    mwfn = find_multiwfn(multiwfn_path)
    _run_multiwfn(mwfn, wfn, _esp_cube_input(spacing), cwd=workdir, timeout=timeout)
    src = os.path.join(workdir, 'mapfunc.cub')
    if not os.path.exists(src):
        raise RuntimeError(f'Multiwfn 未生成 mapfunc.cub（ESP cube 失败）: {src}')
    shutil.copy(src, out_cub)
    return out_cub


def density_cube(wfn: str, out_cub: str, workdir: str = '.',
                 multiwfn_path: str = '', spacing: float = 0.4,
                 timeout: int | None = 3600) -> str:
    """iso 法电子密度网格 cube（与 esp_cube 同一表面格点空间，可直接叠加渲染）。

    spacing: 表面格点间距 Å（默认 0.4，与 esp_cube 保持一致以便叠加）。
    timeout: 单次 Multiwfn 调用超时秒；None=不限。
    返回输出路径。
    """
    mwfn = find_multiwfn(multiwfn_path)
    _run_multiwfn(mwfn, wfn, _esp_density_cube_input(spacing), cwd=workdir,
                  timeout=timeout)
    src = os.path.join(workdir, 'mapfunc.cub')
    if not os.path.exists(src):
        raise RuntimeError(f'Multiwfn 未生成 mapfunc.cub（density cube 失败）: {src}')
    shutil.copy(src, out_cub)
    return out_cub
