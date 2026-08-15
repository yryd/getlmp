#!/usr/bin/env python
"""getlmp 二期 4-site 水冒烟测试（LAMMPS 隐式 M 方案）。

直接运行（无需 pytest）：
    python tests/test_water_4site.py

覆盖：
  T10 tip4pew  4-site 水 + 甲烷（导出层 EP 合并/丢弃、键角修复、tip4p 验证模板）
  T11 tip4p    经典 TIP4P（getlmp 内置 leaprc 分支）
  T12 tip4pd   4-site + 离子（Cl-，验证离子配水链路）
  T13 opc      OPC（最优 4 点水）

验证点：退出码、data.lmp 原子数（EP 丢弃后 = 3×水 + 溶质）、无 EP 原子、
总电荷≈0、键/角数（每水 2 键 1 角）、报告含 4-site 信息与 tip4p pair_style。
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

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


def _cfg_water(name: str, model: str, extra: str = "") -> str:
    """带 1 个甲烷溶质 + 20 个水的 yaml（甲烷 5 原子 + 水 60 原子 = 65）。"""
    return GAFF2_BCC_HEAD.format(name=name) + f"""\
  - smiles: C
    name: methane
    count: 1
    resname: MTH
water:
  model: {model}
  count: 20
packmol:
  box: [0, 0, 0, 20, 20, 20]
{extra}output: data.lmp
workdir: work
seed: 2026
"""


def _parse_header_int(data_path: str, suffix: str) -> int:
    with open(data_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.endswith(suffix):
                return int(line.split()[0])
    raise AssertionError(f"data.lmp 头部未找到 {suffix}: {data_path}")


def _parse_charges_and_types(data_path: str) -> tuple[list[float], list[str]]:
    """读 Atoms 段：返回 (电荷列表, 原子名列表)。"""
    charges: list[float] = []
    names: list[str] = []
    with open(data_path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    i = lines.index("Atoms  # full")
    for line in lines[i + 2:]:
        if not line.strip():
            break
        parts = line.split()
        charges.append(float(parts[3]))
        names.append(parts[8])   # 行尾注释 "# name resname"
    return charges, names


def run_water_test(name: str, model: str,
                   expect_om_qdist: str, extra: str = "",
                   expected_atoms: int = 65) -> None:
    with tempfile.TemporaryDirectory(prefix="getlmp_water_") as tmp:
        yaml_path = os.path.join(tmp, "input.yaml")
        with open(yaml_path, "w", encoding="utf-8") as f:
            f.write(_cfg_water(name, model, extra))

        proc = subprocess.run(
            [sys.executable, str(MAIN_PY), "input.yaml"],
            cwd=tmp, capture_output=True, text=True, timeout=900,
        )
        out = proc.stdout + proc.stderr
        assert proc.returncode == 0, f"getlmp 退出码 {proc.returncode}:\n{out}"

        data_path = os.path.join(tmp, "work", "data.lmp")
        report_path = os.path.join(tmp, "work", "_others", "check_report.txt")
        assert os.path.exists(data_path), f"data.lmp 未生成:\n{out}"
        assert os.path.exists(report_path), f"check_report.txt 未生成:\n{out}"

        # 原子数：甲烷 5 + 20×3 = 65（EPW 已丢弃；含离子时 +离子数）
        atoms = _parse_header_int(data_path, "atoms")
        assert atoms == expected_atoms, \
            f"原子数不符: 期望 {expected_atoms}, 实际 {atoms}\n{out}"

        # 无 EP 原子 / 无 EP 类型
        charges, names = _parse_charges_and_types(data_path)
        assert len(charges) == expected_atoms, \
            f"Atoms 段行数 {len(charges)} != {expected_atoms}\n{out}"
        assert "EPW" not in names, f"data 仍含 EPW 原子: {names}\n{out}"

        # 总电荷 ≈ 0
        total = sum(charges)
        assert abs(total) < 1e-4, f"总电荷 {total:.6f} 不守恒\n{out}"

        # 键/角：甲烷 4 键 6 角 + 水 20×(2 键 1 角) = 44 键 26 角
        nbond = _parse_header_int(data_path, "bonds")
        nangle = _parse_header_int(data_path, "angles")
        assert nbond == 44, f"键数 {nbond} != 44（每水 2 键，EP 键已删）\n{out}"
        assert nangle == 26, f"角数 {nangle} != 26（每水 1 角已补）\n{out}"

        # 报告：4-site 信息 + in 示例引用（4-site 报告不再内联模板，
        # 改为指向 examples/water_4site_lammps_test/<model>/in.<model>.lmp）
        with open(report_path, encoding="utf-8") as f:
            report = f.read()
        assert "结论: 通过" in report, f"校验未通过:\n{report}"
        assert "4-site 虚拟位点: 丢弃 EPW 20 个" in report, \
            f"报告缺 EP 丢弃信息:\n{report}"
        in_rel = (f"examples/water_4site_lammps_test/{model}/"
                  f"in.{model}.lmp")
        assert in_rel in report, \
            f"报告缺 in 示例引用（{in_rel}）:\n{report}"

        # in 示例文件存在且 pair_style 语法正确（类型号按示例体系写死：
        # O=3, H=4, O-H 键 type=2, H-O-H 角 type=2）
        in_path = os.path.join(PROJECT_ROOT, in_rel)
        assert os.path.exists(in_path), f"in 示例缺失: {in_path}"
        with open(in_path, encoding="utf-8") as f:
            in_text = f.read()
        assert (f"pair_style lj/cut/tip4p/long 3 4 2 2 "
                f"{expect_om_qdist} 10.0") in in_text, \
            f"in 示例 pair_style 错误（qdist={expect_om_qdist}）:\n{in_text}"
        assert "kspace_style pppm/tip4p 1e-4" in in_text, \
            f"in 示例 kspace 应 pppm/tip4p:\n{in_text}"


def test_tip4pew() -> None:
    """T10: TIP4P-Ew 4-site 水（标准 AmberTools leaprc 分支）。"""
    run_water_test("water_tip4pew", "tip4pew", "0.1250")


def test_tip4p() -> None:
    """T11: 经典 TIP4P（AmberTools 无 leaprc，走 getlmp 内置分支）。"""
    run_water_test("water_tip4p", "tip4p", "0.1500")


def test_tip4pd_ion() -> None:
    """T12: TIP4P-D + Cl- 离子（离子配水链路 + 4-site 共存）。"""
    # 带 2 个 Cl- 中和 2 个 Na+？—— 用 Na+/Cl- 各 2 保持净电荷 0
    extra = """\
ions:
  - name: Na+
    count: 2
  - name: Cl-
    count: 2
"""
    # 甲烷 5 + 水 60 + 离子 4 = 69；键/角不变（离子无键角）
    run_water_test("water_tip4pd_ion", "tip4pd", "0.1546", extra,
                   expected_atoms=69)


def test_opc() -> None:
    """T13: OPC 4-site 水（电荷分布与 TIP4P 系列不同，仅断言守恒）。"""
    run_water_test("water_opc", "opc", "0.1594")


if __name__ == "__main__":
    tests = [
        ("T10 tip4pew", test_tip4pew),
        ("T11 tip4p", test_tip4p),
        ("T12 tip4pd+离子", test_tip4pd_ion),
        ("T13 opc", test_opc),
    ]
    failed = 0
    for label, fn in tests:
        try:
            fn()
            print(f"  通过 {label}")
        except Exception as e:
            failed += 1
            print(f"  失败 {label}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} 通过")
    sys.exit(1 if failed else 0)
