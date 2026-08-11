#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""getlmp 环境自检：探测 Python 库、外部命令与约定环境变量，输出状态表。

用法：
    python check_env.py              # 全量自检
    python check_env.py --quick      # 只看 QUICK 相关（qm.engine: quick 场景）
    python check_env.py --gaussian   # 只看 Gaussian 相关（qm.engine: gaussian 场景）

退出码：0 = 所有必需项就绪；1 = 存在缺失/无效项（用于 CI 或脚本判断）。

约定（详见 docs/install.md）：
- 所有外部命令（antechamber/parmchk2/tleap/sqm/resp/respgen/packmol/quick/
  g16/formchk/Multiwfn_noGUI）必须在 PATH 中，代码用 shutil.which 探测；
- QUICK 场景需要环境变量 QUICK_BASIS 指向含 6-31G.BAS 的目录；
- Gaussian 场景需要 PATH 含 g16/formchk，且 GAUSS_* 环境变量已配置（.bashrc）。
"""
from __future__ import annotations

import importlib
import os
import shutil
import sys

# (名称, 说明, 必需?) —— 库
PYLIBS = [
    ('rdkit',   'RDKit（SMILES→3D）',          True),
    ('parmed',  'ParmEd（prmtop→data.lmp）',    True),
    ('yaml',    'PyYAML（读 input.yaml）',      True),
    ('numpy',   'NumPy',                        True),
]

# (命令, 说明, 必需?)
# antechamber/parmchk2 链（bcc/abcg2 必需）；tleap/packmol（装盒）；resp/respgen
# （仅 QUICK 旧路径 RESP 需要）；g16/formchk（仅 Gaussian 引擎）；quick（仅 QUICK 引擎）。
CMDS = [
    ('antechamber', 'antechamber（类型+电荷）', True),
    ('parmchk2',    'parmchk2（frcmod 补齐）',  True),
    ('tleap',       'tleap（拓扑+坐标）',       True),
    ('packmol',     'packmol（装盒，count>1）', False),
    ('sqm',         'sqm（AM1-BCC 电荷）',      True),
    ('resp',        'resp（QUICK RESP 拟合）',  False),
    ('respgen',     'respgen（RESP 输入）',     False),
    ('quick',       'quick（QUICK 引擎）',      False),
    ('g16',         'g16（Gaussian 引擎）',     False),
    ('formchk',     'formchk（Gaussian 引擎）', False),
    ('Multiwfn_noGUI', 'Multiwfn（RESP/ESP）',  False),
]

# 环境变量 (名称, 说明, 适用场景, 校验函数)
ENVVARS = [
    ('QUICK_BASIS', 'QUICK 基组数据目录（含 6-31G.BAS）', 'quick',
     lambda v: os.path.isdir(v) and os.path.exists(os.path.join(v, '6-31G.BAS'))),
]


def _check_lib(name: str, desc: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def main() -> int:
    args = sys.argv[1:]
    scope = None
    if '--quick' in args:
        scope = 'quick'
    elif '--gaussian' in args:
        scope = 'gaussian'

    rows: list[tuple[str, str, bool, str]] = []
    for name, desc, required in PYLIBS:
        ok = _check_lib(name, desc)
        rows.append((f'lib  {name}', desc, ok, '' if ok else
                     'pip install rdkit parmed pyyaml numpy（或按 docs/install.md 重建 conda 环境）'))

    for cmd, desc, required in CMDS:
        p = shutil.which(cmd)
        rows.append((f'bin  {cmd}', desc, p is not None,
                     '' if p else '未在 PATH 中，安装与配置见 docs/install.md'))

    for name, desc, when, check in ENVVARS:
        if scope and when != scope:
            continue
        v = os.environ.get(name, '')
        ok = bool(v) and check(v)
        rows.append((f'env  {name}', desc, ok,
                     '' if ok else f'未设置或无效，见 docs/install.md（export {name}=...）'))

    # 打印状态表
    print('=' * 72)
    print('getlmp 环境自检', '(QUICK 场景)' if scope == 'quick'
          else '(Gaussian 场景)' if scope == 'gaussian' else '（全量）')
    print('=' * 72)
    for tag, desc, ok, hint in rows:
        mark = '✓' if ok else '✗'
        print(f'  {mark} [{tag}] {desc}')
        if not ok and hint:
            print(f'      → {hint}')
    print('-' * 72)

    missing = [r for r in rows if not r[2]]
    if not missing:
        print('全部必需项就绪 ✓')
        return 0

    # 区分"必需缺失"与"可选缺失"（可选=对应引擎不使用时不影响）
    required_tags = ({f'bin {c}' for c, _, req in CMDS if req} |
                     {f'lib {n}' for n, _, req in PYLIBS if req})
    required_missing = [r for r in missing if r[0] in required_tags]
    if not required_missing:
        print('提示：存在可选缺失（对应引擎/功能未使用时不影响）。')
        return 0
    print(f'存在 {len(required_missing)} 个必需项缺失，请按 docs/install.md 安装后重跑自检。')
    return 1


if __name__ == '__main__':
    sys.exit(main())
