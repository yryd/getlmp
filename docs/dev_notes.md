# getlmp 开发文档

> 本文档合并原 `docs/` 下四份文档（项目规划、RESP/ReaxFF 验收报告、环境搭建）的精要，
> 作为项目设计与历史经验的总索引。完整原始记录保留在 git 历史（首次提交 c872f75）。

## 1. 项目定位

**一句话**：输入 SMILES 和体系参数，一条命令得到可直接用于 LAMMPS 的 `data.lmp` 与检查报告。

**形态**：`python main.py input.yaml`，不同需求只改 yaml 里的 `forcefield:` / `charge_method:`。

**架构**：分子层（RDKit 3D → antechamber 类型+电荷 → parmchk2）→ 体系层（单分子 tleap /
多分子 packmol 装盒 + tleap loadpdb）→ 导出层（prmtop → data.lmp + 校验 + 报告）。

## 2. 需求矩阵（状态截至 2026-08-10）

| 力场 | 电荷 | 优先级 | 状态 |
|------|------|--------|------|
| GAFF2 | AM1-BCC | 首期 | 已完成（sqm 半经验） |
| ReaxFF | 无（QEq 模拟中算） | 首期 | 已完成（极简 data，集群实跑通过） |
| GAFF2 | RESP | 二期 | 已完成（QUICK 开源 QM，HF/6-31G* 两阶段拟合） |
| GAFF2 | ABCG2 | 二期 | 环境已验证可用（AmberTools 26 内置数据文件） |
| GAFF2 | RESP2 | 三期 | 不支持：antechamber 无 `-c resp2`、QUICK 无隐式溶剂；默认回退 RESP |

## 3. 环境（conda env `smi2data`，2026-08-09 实测）

安装：

```bash
conda create -n smi2data -c conda-forge python=3.12 ambertools=26 rdkit parmed numpy pyyaml -y
```

激活（非交互 shell 必须显式 source）：

```bash
source /home/yryd/packages/miniconda3/etc/profile.d/conda.sh
conda activate /home/yryd/packages/miniconda3/envs/smi2data
```

关键版本：ambertools 26.0（含 antechamber/parmchk2/tleap/sqm/quick/packmol 21.0.1）、
rdkit 2026.03.1、parmed 26.0、numpy 2.4.6、python 3.12。

**环境坑**：

1. **不要单独安装 packmol**：AmberTools 26 自带 packmol 并用约束 `packmol ==9999999999`
   禁止外部安装，`conda install ambertools packmol` 会报 UnsatisfiableError。只装 `ambertools=26`。
2. **WSL 后台任务保活**：`wsl -e bash -c "nohup ... &"` 会在 wsl.exe 退出后杀掉子进程。
   长任务用 `start /b wsl -e bash -c "... & wait"` 保持会话存活。
3. **网络断线恢复**：安装中断报 `CondaHTTPError` 后重跑同一条 `conda create` 命令即可，
   已下载包保留在缓存中，从断点继续。

## 4. RESP 电荷链路（乙醇 + 苯双案例验收通过）

链路：`mol2 → QUICK(HF/6-31G* ESP on vdW 表面) → 自建 .esp → resp 两阶段拟合 → 写回 mol2`。

**三个坑**（均为实测踩过，代码注释已同步）：

1. **QUICK 输入原子行只写 `元素 x y z`（4 字段）**：多写原子序数会被 QUICK 按 4 字段错位读
   （Z 当 x、x 当 y），分子几何错乱、vdW 表面点异常。修复后表面点距原子 1.68–6.40 Å（正常）。
2. **`.esp` 点行顺序 = `ESP X Y Z`**：ESP 值在前、坐标（Bohr）在后，与原子坐标行（仅 3 列）
   区分。写反会导致 resp 拟合发散（电荷量级 10^0~10^2）。
3. **不要用 `antechamber -fi quick` 做等价性**：QUICK 输出无键信息，连接性缺失会导致
   乙醇 9 原子全等效（O 被压到 -0.13、H 全 +0.38，违背化学）。改走 `mol2 → .ac → respgen`，
   苯 6C/6H 等效正确（D6h），乙醇保守退化无假等效。

**验证数据**：乙醇 O=-0.5957、OH-H=+0.3855、CH3/CH2 的 H 各自等效，sum=0.000000，
ESP relative RMS=0.203；苯 6C 全 -0.1252、6H 全 +0.1252（D6h 对称保持），RMS=0.140。
两者均与文献 RESP 电荷一致。

**RESP2 回退方案**：默认回退 RESP（质量好，满足验收）。若博士需要 RESP2（溶剂化 ESP），
需带 PCM 的开源 QM（PSI4/PySCF/ORCA）算水中 HF/6-31G* ESP，新增 ESP 源读入接口即可，
resp 两阶段链路不变。

## 5. ReaxFF 分支（本机 + 集群验收通过）

极简 data（Masses + Atoms，无键项），电荷留 0（QEq 模拟中算），多分子 packmol 装盒。

**集群实测三个坑**（LAMMPS 22 Jul 2025，作业 3251915）：

1. **`atom_style charge` 的 Atoms 段必须 6 列**（`id type q x y z`，无 mol-id）；
   需要分子分组时配置 `reax_atom_style: full`（7 列带 mol-id）。
2. **data 无 Pair Coeffs 段时，`pair_style`/`pair_coeff` 必须放在 `read_data` 之后**
   （与 GAFF2 的"read_data 放所有 style 后"规则互补：data 含 Coeffs 段时反过来）。
3. **`fix qeq/reaxff` 的 params 可用关键字 `reaxff`**：直接从 ffield 提取 QEq 参数，
   无需 param.qeq 文件。

**验证**：乙醇×50 / 40³ 盒 = 450 atoms / 3 types，`read_data` + `run 0` 通过
（VERIFY_EXIT=0），能量与手工修版逐位一致。

## 6. ABCG2（环境已验证，可直接用）

AmberTools 26 已内置 `BCCPARM_ABCG2.DAT` / `ATOMTYPE_ABCG2.DEF`（conda 环境
`$CONDA_PREFIX/dat/antechamber/` 下），`antechamber -c abcg2 -at gaff2` 直接可用，
无需补拷文件。乙醇验证：GAFF2 类型正确（c3/oh/hc/h1/ho），电荷守恒（sum=0.000001）。

## 7. 后续方向（未实现）

- **packmol slab / interface 预设**：slab（双板夹层，GO 夹层接枝）与 interface
  （两相界面，TTSBI-三聚氰氯）的 inp 生成规则、用户可改/不可改字段的设计见
  git 历史中的 `smi2data_plan.md`（3.1 节）。
- **RESP2**：见第 4 节回退方案。

## 8. 交付约定

每个阶段交付：代码 + 检查报告。运行检查报告随每次运行生成于 `workdir/_others/check_report.txt`
；开发经验与验收结论归档至本文档。
