# smi2data 阶段 4 检查报告 — ReaxFF 分支

- 日期: 2026-08-09
- 开发: MDToolAgent ｜ 集群验证: OpsAgent（作业 3251915）
- 状态: ✅ 闭环（本机生成 + 集群 LAMMPS 实跑通过）

## 1. 目标

`forcefield: reaxff` 分支：SMILES → RDKit 3D → 坐标 → 极简 data.lmp
（Masses + Atoms，无键项），多分子 packmol 装盒，电荷留 0（QEq 模拟中算）。

## 2. 实现（模块化扩展）

| 模块 | 改动 |
|------|------|
| `config.py` | `SUPPORTED_FF` 加 `reaxff`；charge_method 须 `none`（QEq 模拟中算）；新增 `reax_elements`（元素顺序→类型号，须与 ffield.reax 一致，默认 C H O N S P F Cl Br I）、`reax_atom_style`（charge/full） |
| `molecule_layer.py` | ReaxFF 分支：跳过 antechamber/parmchk2，SDF → 元素 + 坐标（`_sdf_elements_coords`） |
| `packmol_layer.py` | 新增 `sdf_to_xyz`（ReaxFF 无 mol2，SDF 直接转 packmol 输入） |
| `export_reaxff.py`（新） | 极简 data：Masses（元素质量）+ Atoms（电荷 0）；box 可由 packmol 指定或坐标自动推算 |
| `pipeline.py` | `_run_pipeline_reaxff`：分子层 → packmol → 导出 → 校验（原子数守恒/段计数/类型数）；`REAXFF_IN_TEMPLATE` 验证脚本 |
| `check_lammps_data.py` | parse_data 支持 Atoms 段 6 列（charge）与 7 列（full）双格式自适应 |

## 3. 集群实测暴露并修复的坑（重要）

1. **Atoms 段多输出 molecule-ID 列（7 列）**：
   `atom_style charge` 是 ATOMIC 型，Atoms 段必须 6 列 `id type q x y z`（无 mol-id）。
   修复：默认 charge 6 列；需要分子分组时配置 `reax_atom_style: full`（7 列带 mol-id）。
2. **pair_coeff 必须在 read_data 之后**：
   data 无 Pair Coeffs 段时，`pair_style` + `pair_coeff` 必须放在 `read_data` 之后
   （与阶段 2 的"read_data 放所有 style 后"规则互补：data 含 Coeffs 段时反过来）。
3. **fix qeq/reaxff 的 params 可用关键字 `reaxff`**：直接从 ffield 提取 QEq 参数，
   无需 param.qeq 文件（LAMMPS 22 Jul 2025 实测）。

## 4. 验证结果

### 4.1 本机（conda env smi2data）
- 乙醇×50 / 40³ 盒：450 atoms / 3 atom types，校验全通过
- 乙醇×1：9 atoms，盒自动推算，校验通过
- 回归：GAFF2 苯、乙醇 RESP 均通过（parse_data 双格式自适应后无破坏）

### 4.2 集群（LAMMPS 22 Jul 2025 - Update 4，作业 3251915）
```
Reading data file ...
  reading atoms ...
  450 atoms
  read_data CPU = 0.003 seconds
Reading potential file ffield.reax with DATE: 2011-02-18
   Step    Temp    E_pair      E_mol     TotEng      Press
     0      0     -44369.253   0        -44369.253   7206.8232
VERIFY_EXIT=0
```
- 判据达成：READ_DATA_OK + run 0 无 ERROR
- 能量与手动修版（作业 3251914）逐位一致 → 转换器输出等价

## 5. 案例与配套文件

- `examples/ethanol_reax/`：多分子 bulk（50 乙醇，40³ 盒）
- `examples/ethanol_reax_single/`：单分子（自动盒）
- `examples/ethanol_reax/ffield.reax`：Chenoweth CHO 力场（LAMMPS potentials 官方，
  ffield.reax.cho，元素顺序 C H O）
- `examples/ethanol_reax/in.reax.verify`：集群验证脚本

## 6. 使用方式

```yaml
name: ethanol_reax
forcefield: reaxff
charge_method: none
reax_elements: [C, H, O]        # 须与 ffield.reax 元素顺序一致
molecules:
  - smiles: CCO
    name: ethanol
    count: 50
packmol:
  box: [0, 0, 0, 40, 40, 40]
```

```bash
smi2data input.yaml → data.lmp + check_report.txt
# 集群跑：lmp -in in.reax.verify（read_data → pair_style reaxff → pair_coeff → fix qeq/reaxff）
```

## 7. 已知限制 / 后续

- `reax_elements` 顺序由用户保证与 ffield.reax 一致（工具只做合法性检查：元素必须在列表中）
- 默认 `reax_atom_style: charge` 无分子 ID；需分子分组（如 fix reaxff/bonds 输出）时用 full
- ReaxFF 原子类型更细的映射（同元素多杂化态）留给用户按 ffield.reax 自行扩展
