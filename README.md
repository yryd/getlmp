# getlmp

<!-- 徽章行 -->
[![version](https://img.shields.io/badge/version-0.1.0-blue)](https://github.com/)
[![forcefield](https://img.shields.io/badge/forcefield-GAFF2%20%7C%20ReaxFF-orange)]()
[![charge](https://img.shields.io/badge/charge-AM1--BCC%20%7C%20RESP%20%7C%20RESP2%20%7C%20ABCG2-green)]()
[![status](https://img.shields.io/badge/status-阶段0--5%20闭环-yellow)]()

**一条命令，从 SMILES 到可运行的 LAMMPS `data.lmp`。**

getlmp 是一个分子动力学建模工具链：输入分子的 SMILES 和体系参数（力场、电荷、装盒），
自动完成 3D 建模 → 力场参数 → 电荷拟合 → 装盒 → 拓扑 → LAMMPS data 文件 的全流程，
并输出一份检查报告。不用再手工拼 antechamber / packmol / tleap / ParmEd 的每一步。

```bash
python main.py input.yaml    # → data.lmp + work/check_report.txt
```

---

## Features

- **多力场**：GAFF2/GAFF（经典全原子，GAFF2 参数取自 AmberTools 26 内置 gaff2.dat v2.2.30/2025）与 ReaxFF（反应力场，极简坐标 data）
- **多电荷方案**：AM1-BCC（默认）、RESP / RESP2（Gaussian 16 + Multiwfn 拟合，
  531 路线）、ABCG2；RESP 亦可回退 QUICK（开源）引擎
- **单分子与多分子体系**：单分子自动建盒；多分子用 packmol 装盒，支持任意分子组成
- **完整拓扑**：键/角/二面角/improper 参数由 antechamber + parmchk2 补齐，导出 atom_style full
- **内置校验**：电荷守恒、原子数守恒、段计数、几何自检，跑完即出报告
- **ESP 可视化导出**：单分子 RESP/RESP2 自动产出 iso/pt 两套产物
  （`_others/electrostatic_potential/iso/{density.cub,esp.cub}` 与
  `pt/{mol.pdb,vtx.pdb}`），统一 Multiwfn 导出（Gaussian .fch / QUICK .molden），
  VMD 直接渲染：iso 法密度等值面+ESP 着色，pt 法 Beta 着色
- **集群实测过**：生成的 data.lmp 已多次通过集群 LAMMPS `read_data` + `run 0` 实跑验收
- **可配置**：力场/电荷/装盒参数全部走 yaml，改两行配置即可切换需求

---

## Installation

前置：WSL 或 Linux，已装 Miniconda（getlmp 运行环境用 conda 管理，环境名沿用 smi2data）。

```bash
# 1. 创建 conda 环境（AmberTools 自带 antechamber/parmchk2/tleap/packmol/QUICK）
conda create -n smi2data -c conda-forge python=3.12 \
    ambertools=26 rdkit parmed numpy pyyaml -y
conda activate smi2data

# 2. 运行（无需安装，项目根目录执行）
python main.py input.yaml
```

> 已有环境（本项目当前机器）直接激活：
> `source /home/yryd/packages/miniconda3/etc/profile.d/conda.sh && conda activate smi2data`

---

## Quick Start

准备一个最小配置文件：

```yaml
# input.yaml
name: ethanol
forcefield: gaff2        # 力场：gaff2 / gaff / reaxff
charge_method: bcc       # 电荷：bcc / resp / abcg2（reaxff 时为 none）
molecules:
  - smiles: CCO          # 乙醇
    name: ethanol
    count: 1
output: data.lmp         # 输出 data 文件名（相对 workdir；默认 data.lmp）
workdir: work            # 工作目录
organize_output: true    # 跑完整理：workdir 根目录只留 output 的 lmp + 各分子 mol2 + system.xyz，其余进 _others/
```

运行（项目根目录执行）：

```bash
python main.py input.yaml
```

可选参数：
```bash
python main.py input.yaml --buffer 5.0   # 单分子盒 padding（Å，默认 3.8）
```

每次运行默认导出 `work/system.xyz`（标准 xyz 格式：第 1 行原子数、第 2 行注释、
随后每行 `元素 x y z`），**原子顺序与 data.lmp 完全一致**，可直接用于
OVITO/VMD 可视化或传给其他工具（如 Multiwfn）。

输出（节选）：

```
==== getlmp pipeline: ethanol (gaff2 / bcc, net=0, 单分子) ====
== 分子层: ethanol (CCO) ==
  SMILES → SDF: 9 atoms (incl. H)
  antechamber → mol2 (bcc / gaff2, resname=ETH)
  电荷归一化 → 分子净电荷 = 0
  parmchk2 → frcmod（空：参数全部来自 GAFF 标准库）
== 导出层: prmtop → data.lmp ==
  data.lmp: 9 atoms / 8 bonds / 13 angles / 12 dihedrals / 0 impropers / charge=0.000000

==== 校验结果: 通过 ====
  电荷守恒: sum=0.000000 净电荷=0 通过
  原子数守恒: 9 == 分子层 9 通过
  段计数: 头部与文件记录一致 通过
```

产物：

- `data.lmp` — LAMMPS 可直接 `read_data` 的 data 文件（文件名 = yaml `output:`，默认 `data.lmp`）
- `work/` — 工作目录；`organize_output: true`（默认）时根目录只保留
  `output` 的 lmp + 各分子 `{name}.mol2` + `system.xyz`，其余中间文件（检查报告、prmtop、
  sdf、frcmod、esp.cub、molden、ANTECHAMBER_* 等）自动移入 `work/_others/`；
  同名文件默认直接覆盖只留最新一份，如需保留历史多份可设 `organize_backup: true`

---

## Usage

### 场景 1：多分子体系（苯 100 + 水 1000，bulk 装盒）

```yaml
name: benzene_water
forcefield: gaff2
charge_method: bcc
molecules:
  - smiles: c1ccccc1
    name: benzene
    count: 100
    resname: BEN
  - smiles: O
    name: water
    count: 1000
    resname: WAT
packmol:
  box: [0, 0, 0, 60, 60, 60]    # [xlo, ylo, zlo, xhi, yhi, zhi] Å
```

输出 4200 原子、完整键合拓扑（水 O-H 键齐全），集群 LAMMPS 实跑通过。

### 场景 2：RESP / RESP2 电荷（Gaussian 16 + Multiwfn，531 路线）

```yaml
name: benzene_resp
forcefield: gaff2
charge_method: resp        # resp / resp2（RESP2 需 qm.resp2 + qm.solvent）
net_charge: 0
qm:
  engine: gaussian         # gaussian（默认，G16+Multiwfn）/ quick（QUICK 旧路径）
  method: b3lyp            # 531 推荐 B3LYP-D3(BJ)
  basis: def2TZVP          # G16 命名（def2-TZVP 连字符写法自动兼容）
  # solvent: water         # 空=气相；RESP2 需要溶剂（PCM 单点）
  # resp2: true            # RESP2 开关（需 solvent 非空）
  # delta: 0.5             # RESP2 δ 混合系数
  # multiwfn_path: ''      # 空=自动探测（~/packages/soft/multiwfn/）
molecules:
  - smiles: c1ccccc1
    name: benzene
    count: 1
# esp 段可调 ESP 导出（默认开启，仅单分子 RESP/RESP2 生效）：
# esp:
#   enabled: true   # 不需要 ESP 可视化时改为 false
#   # spacing/buffer 已废弃（旧点电荷近似 cube 参数），Multiwfn 网格自动
```

流程：SMILES → 3D → antechamber（bcc 类型占位）→ g16 单点（pop=MK ESP）
→ formchk → Multiwfn RESP 拟合 → 电荷写回 mol2 → tleap → data.lmp。

单分子 RESP/RESP2 时自动多产出可视化文件（`_others/electrostatic_potential/` 内）：
- `iso/density.cub` + `iso/esp.cub`：iso 法——Multiwfn 严格 QM 电子密度与 ESP 网格
  （同一网格，VMD 用 density 等值面 0.001 画表面 + esp 着色）
- `pt/mol.pdb` + `pt/vtx.pdb`：pt 法——分子结构与密度 0.001 闭合等值面顶点
  （B 因子字段存 ESP kcal/mol，VMD Beta 着色直接看图）
- `{name}.fch` / `{name}.molden`：波函数归档（Gaussian / QUICK，二次分析用）

RESP 亦可回退 QUICK 引擎（`qm.engine: quick`，HF/6-31G* + resp 两阶段拟合，
Multiwfn 读 molden 产出同样的 iso/pt 两套产物），示例见 `examples/PIP.yaml`。

### 场景 3：ReaxFF 反应力场（坐标 + 元素，无键项，QEq 模拟中算）

```yaml
name: ethanol_reax
forcefield: reaxff
charge_method: none
reax_elements: [C, H, O]   # 元素顺序 = data 原子类型号，须与 ffield.reax 一致
molecules:
  - smiles: CCO
    name: ethanol
    count: 50
packmol:
  box: [0, 0, 0, 40, 40, 40]
```

输出极简 data（Masses + Atoms），配 `ffield.reax` + `pair_style reaxff` 即可运行，集群实测通过。

> 可直接运行的完整示例见 `examples/`：`cd examples && python ../main.py methane.yaml`；
> 注意：产物生成在 `examples/work/`（临时目录，已在 .gitignore，跑完可删）；
> 设计与经验总结见 `docs/dev_notes.md`。

---

## 输出与校验

每次运行自动完成 5 项校验，全部通过才报告：

| 校验 | 说明 |
|------|------|
| 电荷守恒 | Σq ≈ 净电荷（容差 0.01） |
| 原子数守恒 | 输入拷贝数 × 每分子原子数 == data 原子数 |
| 段计数 | 头部声明 vs 文件实际记录数 |
| 类型计数 | 头部 vs Coeffs 段 |
| 拓扑期望 | 多分子键/角数 == Σ count × 单分子拓扑（键合完整性） |

生成的 data.lmp 同时经**集群 LAMMPS 实跑验收**（`read_data` + `run 0`，VERIFY_EXIT=0）——
验收记录见 `docs/dev_notes.md`。

### 流程会产出哪些文件（workdir 内）

> `organize_output: true`（默认）跑完自动整理：根目录只留 `output` 的 lmp + 各 `{name}.mol2` + `system.xyz`，
> 下面表格里的其余文件全部移入 `work/_others/`。

| 类别 | 文件 | 说明 |
|------|------|------|
| **主产物** | `data.lmp` | 最终 LAMMPS data（交付物） |
| | `system.xyz` | 体系标准 xyz（原子顺序与 data.lmp 一致，OVITO/VMD/Multiwfn 用） |
| | `check_report.txt` | 检查报告（交付物） |
| 中间必要 | `{name}.sdf` | RDKit 3D 构象（SMILES→3D） |
| | `{name}.mol2` | antechamber 产物：GAFF 类型+坐标+键+电荷（**任何分支默认都有**） |
| | `{name}.frcmod` | parmchk2 补参（空 = GAFF 标准库齐全） |
| | `{name}.prmtop` / `{name}.inpcrd` | tleap 拓扑 + 坐标（→ data.lmp 的原料） |
| | `tleap.in` / `leap.log` | tleap 输入/日志 |
| | `sqm.in/out/pdb` | AM1-BCC 半经验计算（仅 BCC/ABCG2） |
| RESP 额外 | `{name}_quick.in/.out/.vdw` | QUICK 输入/日志/vdW 表面 ESP 点 |
| | `_resp_tmp_{name}/` | resp 两阶段拟合中间目录 |
| | `_others/electrostatic_potential/iso/` | iso 法：density.cub + esp.cub（Multiwfn 严格 QM） |
| | `_others/electrostatic_potential/pt/` | pt 法：mol.pdb + vtx.pdb（B 因子=ESP） |
| | `{name}.fch` / `{name}.molden` | 波函数（可视化/二次分析） |
| 噪音残留 | `ANTECHAMBER_*.AC`、`ATOMTYPE.INF` | antechamber 临时文件，可忽略/删除 |

---

## Documentation

| 文档 | 内容 |
|------|------|
| [`examples/`](examples/) | 可直接运行的示例配置（甲烷/乙醇单体 + 混合体系 + RESP/RESP2 + QUICK 回退） |
| [`docs/yaml_config.md`](docs/yaml_config.md) | **YAML 配置参考**：全部选项、类型、默认值、校验规则、组合约束 |
| [`docs/dev_notes.md`](docs/dev_notes.md) | 开发文档：规划与需求矩阵、环境搭建、技术坑、验收结论、后续方向 |

---

## FAQ

**Q: 支持哪些力场和电荷？**
GAFF2/GAFF + AM1-BCC（默认）、RESP、RESP2、ABCG2；ReaxFF（QEq 电荷模拟中算）。
RESP/RESP2 的 QM 引擎默认 Gaussian 16（`qm.engine: gaussian`，需 WSL 已装
G16 与 Multiwfn，安装与 PATH 配置见 [`docs/install.md`](docs/install.md)）；
未装时可回退 `qm.engine: quick`（QUICK 免费，但无 RESP2、无隐式溶剂）。
安装好环境后可用 `python check_env.py` 自检。

**Q: 多分子怎么装盒？**
`count > 1` 自动启用 packmol（需在 `packmol.box` 给出盒尺寸）。
当前支持 bulk 预设；slab / interface 在规划中。

**Q: ReaxFF 的原子类型怎么定的？**
data 的原子类型号 = `reax_elements` 列表顺序（默认 C H O N S P F Cl Br I）。
该顺序必须与你使用的 `ffield.reax` 元素顺序一致，工具会做元素存在性检查。

**Q: 需要本地装 LAMMPS 吗？**
不需要。工具只产出 data 文件；LAMMPS 可读性验收在集群完成（本项目已实测 4 次作业）。

**Q: 生成的 data 能直接用吗？**
能。GAFF2 输出为 `atom_style full` 完整 data；ReaxFF 输出为 `atom_style charge`
极简 data（Masses + Atoms）。两者都有配套的 LAMMPS 验证模板（见检查报告）。

---

## Contributing（开发）

```text
main.py                  # 标准入口：python main.py input.yaml
src/
├── config.py            # yaml 配置解析与校验
├── pipeline.py          # 流水线编排（GAFF2 主链 / ReaxFF 分支）+ 检查报告
├── molecule_layer.py    # SMILES → 3D → antechamber → parmchk2；RESP/RESP2（G16+Multiwfn 或 QUICK）
├── packmol_layer.py     # 多分子装盒
├── system_layer.py      # tleap 拓扑 ×N + 坐标 → prmtop
├── export_layer.py      # prmtop → data.lmp（atom_style full）
├── export_reaxff.py     # 坐标/元素 → data.lmp（atom_style charge）
├── multiwfn.py          # Multiwfn RESP/RESP2 拟合 + ESP 可视化导出（iso/pt, fch/molden）
├── check_lammps_data.py # 不依赖 LAMMPS 的格式自检
├── check_env.py         # 环境自检（库/命令/环境变量，见 docs/install.md）
└── prmtop_to_lammps.py  # Amber → LAMMPS 参数转换
```

测试（无需 pytest，直接运行；每个用例在临时目录执行，不污染项目）：

```bash
python tests/test_build_systems.py    # 甲烷/乙醇单体 + 混合体系，3/3
```

改动后回归：运行上述测试，或在临时目录准备 `input.yaml` 执行 `python main.py input.yaml` 验证输出与校验报告。

每个阶段交付：代码 + 检查报告（运行报告在 `workdir/_others/`，开发经验归档 `docs/dev_notes.md`）。

---

## License

内部科研工具（膜分离/水处理 + 计算化学课题组），暂未对外发布。
分子结构数据与项目文件均在保密范围，请勿外传。
