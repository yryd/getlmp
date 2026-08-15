#!/usr/bin/env python
"""getlmp 冒烟测试：甲烷 / 乙醇 单分子 + 混合体系。

直接运行（无需 pytest）：
    python tests/test_build_systems.py

入口用项目内 `main.py`（无需安装）。
每个用例都在临时目录中执行，产物（work/、data.lmp）不会写入项目目录。
验证点：命令退出码、data.lmp 原子数、检查报告校验结论（通过）。
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

# 避免运行测试时在项目 src/ 生成 __pycache__
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MAIN_PY = PROJECT_ROOT / "main.py"

GAFF2_BCC_HEAD = """\
name: {name}
forcefield: gaff2
charge_method: bcc
net_charge: 0
molecules:
"""


def _cfg_single(name: str, smiles: str) -> str:
    return GAFF2_BCC_HEAD.format(name=name) + f"""\
  - smiles: {smiles}
    name: {name}
    count: 1
output: data.lmp
workdir: work
seed: 2026
"""


def _cfg_mix() -> str:
    return GAFF2_BCC_HEAD.format(name="methane_ethanol_mix") + """\
  - smiles: C
    name: methane
    count: 10
    resname: MTH
  - smiles: CCO
    name: ethanol
    count: 5
    resname: ETH
packmol:
  box: [0, 0, 0, 30, 30, 30]
output: data.lmp
workdir: work
seed: 2026
"""


def _parse_atoms(data_path: str) -> int:
    """从 data.lmp 头部解析原子数（'N atoms' 行）。"""
    with open(data_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.endswith("atoms"):
                return int(line.split()[0])
    raise AssertionError(f"data.lmp 头部未找到 atoms 行: {data_path}")


def run_getlmp(yaml_text: str, expected_atoms: int) -> None:
    """在临时目录运行 main.py 入口，断言退出码 / 原子数 / 校验通过。"""
    with tempfile.TemporaryDirectory(prefix="getlmp_test_") as tmp:
        yaml_path = os.path.join(tmp, "input.yaml")
        with open(yaml_path, "w", encoding="utf-8") as f:
            f.write(yaml_text)

        proc = subprocess.run(
            [sys.executable, str(MAIN_PY), "input.yaml"],
            cwd=tmp,
            capture_output=True,
            text=True,
            timeout=600,
        )
        out = proc.stdout + proc.stderr
        assert proc.returncode == 0, f"getlmp 退出码 {proc.returncode}:\n{out}"

        data_path = os.path.join(tmp, "work", "data.lmp")
        xyz_path = os.path.join(tmp, "work", "system.xyz")
        # organize_output 默认开启：data.lmp/mol2/system.xyz 留根目录，其余（含检查报告）进 _others/
        report_path = os.path.join(tmp, "work", "_others", "check_report.txt")
        assert os.path.exists(data_path), f"data.lmp 未生成:\n{out}"
        assert os.path.exists(report_path), f"check_report.txt 未生成（应在 _others/）:\n{out}"
        assert os.path.exists(xyz_path), f"system.xyz 未生成（应默认导出）:\n{out}"
        assert not os.path.exists(os.path.join(tmp, "work", "_others", "system.xyz")), \
            f"system.xyz 不应被移入 _others/:\n{out}"

        atoms = _parse_atoms(data_path)
        assert atoms == expected_atoms, (
            f"原子数不符: 期望 {expected_atoms}, 实际 {atoms}\n{out}"
        )

        # system.xyz 原子数须与 data.lmp 一致（原子顺序一致，可喂 OVITO/VMD/Multiwfn）
        with open(xyz_path, encoding="utf-8") as f:
            xyz_n = int(f.readline().split()[0])
        assert xyz_n == atoms, f"system.xyz 原子数 {xyz_n} != data.lmp {atoms}\n{out}"

        with open(report_path, encoding="utf-8") as f:
            report = f.read()
        assert "结论: 通过" in report, f"校验未通过:\n{report}"


def test_methane_single() -> None:
    """甲烷单体：CH4 = 5 原子。"""
    run_getlmp(_cfg_single("methane", "C"), expected_atoms=5)


def test_ethanol_single() -> None:
    """乙醇单体：C2H6O = 9 原子。"""
    run_getlmp(_cfg_single("ethanol", "CCO"), expected_atoms=9)


def test_methane_ethanol_mix() -> None:
    """混合体系：甲烷 10 + 乙醇 5 → 10*5 + 5*9 = 95 原子（packmol 装盒）。"""
    run_getlmp(_cfg_mix(), expected_atoms=95)


def test_custom_inp_reuse() -> None:
    """自定义 packmol inp（乙醇拆 1 fixed + 4 free）+ reuse_molecule 复用分子层。

    流程：第一次跑生成默认 inp → 改造（fixed 钉在盒中心）→ yaml 加
    inp_file + reuse_molecule → 第二次跑。验证：退出码、[reuse] 复用提示、
    [custom inp] 生效、95 原子守恒、校验通过。
    """
    with tempfile.TemporaryDirectory(prefix="getlmp_test_") as tmp:
        yaml_path = os.path.join(tmp, "input.yaml")
        with open(yaml_path, "w", encoding="utf-8") as f:
            f.write(_cfg_mix())

        # 第一次跑：生成默认 packmol.inp（留在 workdir 根目录）
        r1 = subprocess.run(
            [sys.executable, str(MAIN_PY), "input.yaml"], cwd=tmp,
            capture_output=True, text=True, timeout=600)
        assert r1.returncode == 0, r1.stdout + r1.stderr
        inp_path = os.path.join(tmp, "work", "packmol.inp")
        assert os.path.exists(inp_path), f"默认 packmol.inp 未保留在 workdir:\n{r1.stdout + r1.stderr}"
        with open(inp_path, encoding="utf-8") as f:
            default_inp = f.read()

        # 改造：乙醇块 → 1 fixed（盒中心）+ 4 自由
        # （write_inp 的 inside box 已带 tolerance/2=1.0 Å 缩进，见 packmol_layer.py）
        ethanol_block = ("structure ethanol.xyz\n"
                         "  number 5\n"
                         "  inside box 1.000 1.000 1.000 29.000 29.000 29.000\n"
                         "end structure\n")
        assert ethanol_block in default_inp, "默认 inp 未找到乙醇块"
        fixed_block = ("structure ethanol.xyz\n"
                       "  number 1\n"
                       "  fixed 15.0 15.0 15.0 0.0 0.0 0.0 1.0\n"
                       "end structure\n"
                       "\n"
                       + ethanol_block.replace("number 5", "number 4"))
        my_inp = os.path.join(tmp, "work", "my.inp")
        with open(my_inp, "w", encoding="utf-8") as f:
            f.write(default_inp.replace(ethanol_block, fixed_block))

        # 第二次跑：自定义 inp + 分子层复用
        cfg2 = _cfg_mix().replace(
            "packmol:\n  box: [0, 0, 0, 30, 30, 30]",
            "packmol:\n  box: [0, 0, 0, 30, 30, 30]\n  inp_file: work/my.inp",
        ) + "reuse_molecule: true\n"
        with open(yaml_path, "w", encoding="utf-8") as f:
            f.write(cfg2)
        r2 = subprocess.run(
            [sys.executable, str(MAIN_PY), "input.yaml"], cwd=tmp,
            capture_output=True, text=True, timeout=600)
        out2 = r2.stdout + r2.stderr
        assert r2.returncode == 0, f"第二次跑退出码 {r2.returncode}:\n{out2}"
        assert "[reuse] methane" in out2, f"甲烷未复用分子层:\n{out2}"
        assert "[custom inp]" in out2, f"未使用自定义 inp:\n{out2}"

        data_path = os.path.join(tmp, "work", "data.lmp")
        assert _parse_atoms(data_path) == 95, f"第二次跑原子数不符:\n{out2}"
        report_path = os.path.join(tmp, "work", "_others", "check_report.txt")
        with open(report_path, encoding="utf-8") as f:
            report = f.read()
        assert "结论: 通过" in report, f"第二次跑校验未通过:\n{report}"

        # 自定义 inp 不被 organize 移走（下次还能用）
        assert os.path.exists(my_inp), "自定义 inp 被 organize 移走了"


def main() -> int:
    tests = [
        ("甲烷单体", test_methane_single),
        ("乙醇单体", test_ethanol_single),
        ("混合体系(甲烷10+乙醇5)", test_methane_ethanol_mix),
        ("自定义inp+分子层复用", test_custom_inp_reuse),
    ]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  通过 {name}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"  失败 {name}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} 通过")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
