# smi2data 流水线建设计划
**SMILES → LAMMPS data 一键自动化（WSL 本地）**

作者：SimAgent | 日期：2026-08-08 | 状态：待评审

---

## 1. 项目目标

**一句话**：给出 SMILES 和体系参数，一条命令得到可直接用于 LAMMPS 的 `data.lmp` 文件。

**对用户的形态**（唯一入口，不换工具）：
```
smi2data input.yaml  →  data.lmp + 检查报告
```
不同需求只改配置文件里的 `forcefield:` / `charge:` 两行。

---

## 2. 需求矩阵

| # | 力场 | 电荷 | 分子建模 | 类型指认 | 电荷计算 | 优先级 | 状态 |
|---|------|------|---------|---------|---------|--------|------|
| 1 | GAFF2 | AM1-BCC | RDKit | antechamber | sqm（内置） | ★★★ | 首期交付 |
| 2 | GAFF2 | RESP | RDKit | antechamber | QUICK（开源 QM，AmberTools26 支持） | ★★ | 二期 |
| 3 | ReaxFF | 无（QEq 模拟中算） | RDKit | 元素/ReaxFF 类型 | 不需要 | ★★★ | 首期交付 |
| 4 | GAFF2 | ABCG2 | RDKit | antechamber | `-c abcg2`（Amber 官方源码已集成） | ★★ | 阶段0验证后可用 |
| 5 | GAFF2 | RESP2 | RDKit | antechamber | AmberTools 直接支持（用户调研，阶段3实测） | ★★ | 并入阶段 3 验证 |

> **ABCG2 调研结论（2026-08-08，已实证）**：
> - Amber 官方 GitHub 仓库（Amber-MD/AmberClassic）`src/antechamber/charge.c` **已包含 `-c abcg2` 分支**，配套数据文件 `dat/antechamber/BCCPARM_ABCG2.DAT` 与 `ATOMTYPE_ABCG2.DEF`（ABCG2 独立原子类型定义）均存在
> - 该支持由 Amber 核心开发者 David A. Case 于 **2024-10-21** 提交（晚于 AmberTools24 发布，故 24 官方版不含；AmberTools26 于 2026-04-30 发布，**大概率已包含**）
> - ⚠️ 最终以环境实测为准：阶段 0 安装后直接跑 `antechamber -c abcg2` 验证；若缺失，从 AmberClassic 仓库补拷 `BCCPARM_ABCG2.DAT`/`ATOMTYPE_ABCG2.DEF` 或编译该 antechamber

---

## 3. 架构设计

```
┌─ 用户层 ────────────────────────────────┐
│  input.yaml + `smi2data` 一条命令        │
├─ 流水线编排层（Python，唯一入口）──────────┤
│  解析配置 → 调度各步骤 → 校验 → 报告      │
├─ 分子层 ────────────────────────────────┤
│  RDKit: SMILES → 3D (SDF)               │
│  antechamber: 类型(gaff2) + 电荷(bcc/resp) │
│  parmchk2: 缺失参数补齐 → frcmod          │
├─ 体系层 ────────────────────────────────┤
│  packmol: 装盒 → 坐标 PDB                │
│  tleap: 分子拓扑×N 拷贝 + 坐标绑定 → prmtop │
├─ 导出层 ────────────────────────────────┤
│  ParmEd: prmtop → data.lmp (atom_style full) │
└────────────────────────────────────────┘
```

**MSTK 定位**（后续可选层，不阻塞主流程）：
- 主流程建模/参数走 AmberTools（官方保证）
- MSTK 用于：ommhelper（slab 修正/电场/壁势）、SFE 自由能、scheduler（Slurm 高通量）
- 接入条件：GAFF2 zff 转换器 + mol2 电荷回填脚本（阶段 5）

---

### 3.1 packmol 预设驱动设计（已与用户敲定 2026-08-08）

**原则**：流水线只生成 inp 模板，用户在其基础上手工修改几何控制，改完再跑。**数量与分子增删不做**（number 由流水线按 yaml 生成，用户不改）。

**工作流（三阶段分离）**：
```
smi2data phase1   # 只生成：yaml → 预设 inp + 单分子 xyz 骨架（不跑 packmol）
   ↓  用户手工改 inp（region / fixed / constrain / seed / tolerance）
smi2data phase2   # 只运行：把当前 inp 丢给 packmol → 坐标 _packed.xyz
   ↓
smi2data phase3   # 绑定拓扑：坐标按 structure 顺序 ↔ 分子拷贝（数量以 yaml 为准）→ data.lmp
```

**yaml 结构**：
```yaml
molecules:
  - smiles: O                 # 或 file: xxx.xyz（用户提供的板/大分子）
    count: 2000
  - smiles: CCO
    count: 100
packmol:
  preset: bulk                # bulk / slab / interface
  box: [0, 0, 0, 60, 60, 60]
  seed: 2026
  tolerance: 2.0              # Å
```

**三个预设的 inp 生成规则**：

1. **bulk（均相，默认）** — 每个分子一个 structure 块，`inside box` 整盒
2. **slab（双板夹层，GO 夹层接枝）** — 前两个分子视为上下板：`fixed` 锁位置+朝向（6 值）；其余分子 `inside box` 限中间 z 域 `[z_lo, z_hi]`
   ```yaml
   packmol:
     preset: slab
     box: [0, 0, 0, 50, 50, 20]
     slab_z: [7, 13]          # 中间填充域（Å）
   ```
   → inp：
   ```
   structure go_sheet1.xyz
     number 1
     fixed 25 25 5 0 0 0
   end structure
   structure go_sheet2.xyz
     number 1
     fixed 25 25 15 0 0 0
   end structure
   structure _MO_2.xyz        # O
     number 800
     inside box 0 0 7 50 50 13
   end structure
   structure _MO_3.xyz        # 两性离子
     number 20
     inside box 0 0 7 50 50 13
   end structure
   ```
3. **interface（两相界面，TTSBI-三聚氰氯）** — 分子按 `phase: A/B` 分组，同相共享 z 域（每相可放多种分子），界面在 `interface_z`
   ```yaml
   molecules:
     - smiles: O
       count: 1000
       phase: A
     - smiles: CCO
       count: 100
       phase: A
     - smiles: c1ccccc1
       count: 200
       phase: B
   packmol:
     preset: interface
     box: [0, 0, 0, 60, 60, 60]
     interface_z: 30          # 界面位置（Å）
   ```
   → inp：A 相分子 `inside box 0 0 0 60 60 30`，B 相分子 `inside box 0 0 30 60 60 60`（用户可改 interface_z、加 phase C、加 fixed 板等）

**用户可改 / 不可改**：
- ✅ 可改：region、fixed、constrain_distance/rotation、seed、tolerance、interface_z、slab_z——只影响坐标落点，不影响绑定
- ⚠️ 不改：number、分子增删——拓扑拷贝与坐标绑定以 yaml 为准，改 inp 的 number 会导致 phase3 校验失败

**用户输入文件**：slab 的板（`go_sheet*.xyz`）、interface 中需固定的大分子用 `file:` 指定，流水线不生成。

**与 MSTK 的关系**：MSTK 自带 `scale_with_packmol` 只支持 inside box/slab 分层，无 fixed/区域组合 → 流水线**自己生成 inp**（控制面全开），坐标回填仍可复用 MSTK 的 `update_molecules + set_positions` 思路，导出用 ParmEd/AmberTools 主流程。

---

## 4. 环境准备（WSL2，Ubuntu 26.04 已装）

```
1. wsl --install 状态确认 ✅（Ubuntu-26.04, WSL2, 当前 Stopped）
2. 启动 WSL，更新系统包
3. 安装 Miniconda（Linux x86_64）
4. conda 创建专用环境，安装：
   conda install -c conda-forge ambertools=26 packmol rdkit parmed numpy pandas
   # 说明：
   #  - ambertools: antechamber/parmchk2/tleap/sqm/QUICK（QUICK 用于 RESP）
   #  - packmol: 装盒
   #  - rdkit: SMILES→3D
   #  - parmed: prmtop→LAMMPS data
5. 验证：antechamber --help / parmchk2 --help / packmol / tleap 均可运行
```

---

## 5. 实施阶段（含验收标准）

### 阶段 0：环境搭建 ✅ 目标：工具链全部可运行
- [ ] WSL 启动 + conda 安装 + 工具安装
- [ ] 验收：`antechamber -h`、`parmchk2 -h`、`packmol`、`tleap -h` 无报错；`python -c "import rdkit, parmed"` 成功
- [ ] **ABCG2 实测**：跑 `antechamber -i test.sdf -fi sdf -o test.mol2 -fo mol2 -c abcg2 -at gaff2 -nc 0`，检查是否报 "unknown charge type"；同时确认 `BCCPARM_ABCG2.DAT` 是否随包安装（`find $CONDA_PREFIX -name "BCCPARM*"`）
- [ ] 若 ABCG2 缺失：从 AmberClassic 仓库下载 `BCCPARM_ABCG2.DAT` + `ATOMTYPE_ABCG2.DEF` 放入 antechamber 数据目录（记录处理方式）

### 阶段 1：MVP — 单分子 GAFF2 + AM1-BCC → data.lmp
目标：最小闭环，先跑通再扩展
- [ ] 编写 `smi2data` 骨架（CLI + yaml 解析）
- [ ] 分子层：RDKit SMILES→SDF → antechamber(bcc,gaff2) → mol2 → parmchk2 → frcmod
- [ ] 导出层：单分子 → ParmEd → data.lmp
- [ ] 验收：苯/乙醇各出 data.lmp；**电荷总和=净电荷**；原子数正确；LAMMPS 能读入（pair_style lj/cut/coul/long + 对应 data）
- [ ] 产出：`mstk_examples/` 目录下测试案例

### 阶段 2：多分子体系（packmol + tleap）
- [ ] 体系层：packmol 装盒 → 坐标 PDB
- [ ] tleap：多分子 loadmol2 + 坐标 → 体系 prmtop/inpcrd
- [ ] 导出：ParmEd 转 data.lmp
- [ ] 验收：苯100+水1000 体系；原子数守恒；水分子键合完整；电荷守恒

### 阶段 3：RESP / RESP2 分支
- [ ] antechamber + QUICK 算 RESP 电荷（AmberTools26 特性）
- [ ] 实测 `-c resp2` 是否可用（用户调研：AmberTools 直接支持；若版本不支持则记录回退方案）
- [ ] 验收：与 AM1-BCC 输出结构一致，仅电荷不同；RESP/RESP2 电荷总和=净电荷

### 阶段 4：ReaxFF 分支
- [ ] RDKit 3D → 坐标 → 极简 data（Masses + Atoms，无键项）
- [ ] 多分子装盒（坐标层面）
- [ ] 验收：data 可配 ffield.reax 运行；电荷可留 0（QEq 计算）

### 阶段 5：ABCG2 分支（环境验证通过后与 AM1-BCC 同路径）
- [ ] 配置 `charge: abcg2` → antechamber `-c abcg2`（若阶段 0 确认缺失，先补数据文件/编译）
- [ ] 验证：与 AM1-BCC 输出对比（类型一致、电荷不同）；FreeSolv 类简单分子电荷守恒
- [ ] 产出：ABCG2 与 AM1-BCC 的电荷差异报告（供物理分析参考）

### 阶段 6：MSTK 接入（可选层）
- [ ] 编写 `gaff2.dat → zff` 转换器
- [ ] 编写 mol2 电荷回填脚本（元素+键连图同构匹配，不按顺序硬填）
- [ ] 验证：MSTK System 导出 LAMMPS 与 ParmEd 结果对比（能量一致性）
- [ ] 验收：同一分子 MSTK/ParmEd 两路 data 能量一致（≤1e-3 相对差）

### 阶段 7：~~RESP2 暂缓~~ → 已并入阶段 3（用户调研：AmberTools 直接支持）
- [ ] 归档：RESP2 验证结果与命令行记录写入阶段 3 验收报告

---

## 6. 验证方案（贯穿所有阶段）

| 检查项 | 方法 |
|--------|------|
| 电荷守恒 | sum(charge) ≈ 净电荷（容差 0.01） |
| 原子数守恒 | 输入 counts × 每分子原子数 = data 原子数 |
| 拓扑完整性 | 键/角/二面角数量合理；parmchk2 无 missing（或已有 frcmod 补齐） |
| 能量一致性 | 单点能：两套独立导出结果对比；与 OpenMM 读 prmtop 交叉验证 |
| LAMMPS 可运行 | 最小 in.lmp 能启动并跑 0 步 |

---

## 7. 风险与应对

| 风险 | 影响 | 应对 |
|------|------|------|
| WSL 内 conda-forge 安装冲突 | 中 | 独立 conda env；必要时装系统级依赖 |
| RESP(QUICK) 对新分子收敛差 | 中 | 记录失败案例，回退 AM1-BCC 并提示 |
| GAFF2 缺参数（复杂杂环） | 中 | parmchk2 自动补；报告列出所有补充项 |
| ABCG2 代码拿不到 | 低 | 已预留接口，不影响主线 |
| ParmEd 转 LAMMPS 的 1-4 缩放/非键混合规则 | 高 | 导出后人工抽查 data 的 pair_coeff 与 1-4 缩放，与 prmtop 数值核对 |

---

## 8. 交付物

1. `smi2data` 脚本包（CLI + 模块化步骤 + yaml 配置模板）
2. 测试案例集（苯、乙醇、水、多分子体系）
3. 使用文档（README：配置说明 + 常见问题）
4. 阶段 5 完成后：GAFF2 zff + MSTK 接入示例
5. 每个阶段的检查报告

---

## 9. 里程碑（估计）

| 阶段 | 内容 | 估时 |
|------|------|------|
| 0 | 环境搭建 + ABCG2 实测 | 0.5 天 |
| 1 | MVP 单分子 | 1 天 |
| 2 | 多分子体系 | 1 天 |
| 3 | RESP 分支 | 1 天 |
| 4 | ReaxFF 分支 | 0.5 天 |
| 5 | ABCG2 分支 | 0.5 天 |
| 6 | MSTK 接入 | 2 天 |
| 7 | 交付文档（RESP2 结论并入阶段 3 报告） | — |

**待用户确认**：阶段 4（ReaxFF）的 ffield.reax 力场文件是否已有目标元素的参数？阶段 6（MSTK 接入）是否本期要做，还是先跑通 1-5 再说？
