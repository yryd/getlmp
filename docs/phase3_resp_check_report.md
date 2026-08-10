# 阶段 3 验收检查报告：RESP 电荷链路（QUICK + resp 两阶段）

- 日期：2026-08-09
- 执行：MDToolAgent（smi2data 开发）
- 状态：✅ 通过（乙醇 + 苯 双案例）

## 1. 链路

SMILES → SDF(RDKit) → antechamber(bcc 类型占位) → mol2 → QUICK(HF/6-31G* ESP on vdW 表面)
→ antechamber(mol2→.ac) + respgen(RESP1.IN/RESP2.IN) → 自建 .esp（espgen 兼容格式）
→ resp 两阶段拟合 → 电荷写回 mol2 → 下游（parmchk2/tleap/ParmEd）不变。

## 2. 修复的三个坑（乙醇电荷异常根因）

### 坑 1：QUICK 输入原子行多写原子序数
- 现象：`C 6 0.85 0.37 -0.153`（5 字段）→ QUICK 按 4 字段读，把 Z 当 x、x 当 y…
  分子几何错乱（C1 被放到 x=6.0 Å），vdW 表面点距原子 0.2 Å、X 范围到 11 Å。
- 修复：`元素 x y z`（QUICK READ_COORD 只读 4 字段）。
- 验证：vdW 点距原子 1.68–6.40 Å（正常表面），点数 17112→9777。

### 坑 2：.esp 点行顺序
- 现象：写 `X Y Z Q`（坐标在前）→ resp 拟合发散（电荷量级 10^0~10^2，如 C=-6.7）。
- 真相：espgen/resp 约定点行 = `ESP X Y Z`（ESP 值在前、坐标 Bohr 在后），
  与原子坐标行（仅 3 列）区分。
- 修复：点行改回 `值 坐标` 顺序。

### 坑 3：antechamber `-fi quick` 的等价性 bug
- 现象：respgen 从 -fi quick 的连接性（quick 输出无键信息）判断等价性，
  乙醇 9 原子全等效 → O 电荷被压到 -0.13、H 全部 +0.38（违背化学）。
- 修复：改走 mol2 → .ac（连接性完整）→ respgen。
  - 苯：6C/6H 等效正确（D6h）
  - 乙醇：resp1 保守全独立；resp2 阶段甲基 3H、亚甲基 2H 分别等效（正确）
  - 无假等效（RDKit CanonicalRankAtoms 有假阳性，已弃用该路径）

## 3. 验证结果

### 乙醇（CCO，GAFF2，net=0）
| 原子 | 电荷 |
|---|---|
| C1 (CH3) | +0.0102 |
| C2 (CH2) | +0.0812 |
| O1 | -0.5957 |
| H1-H3 (CH3) | +0.0034（等效 ✓） |
| H4-H5 (CH2) | +0.0543（等效 ✓） |
| H6 (OH) | +0.3855 |
- sum = 0.000000 ✓；ESP relative RMS = 0.203
- 化学合理：O 负、OH-H 正、烷基 H 近 0

### 苯（c1ccccc1，GAFF2，net=0）
- 6 个 C 全部 -0.1252（max diff 0.00000）✓ D6h 对称保持
- 6 个 H 全部 +0.1252 ✓
- sum = 0.000000 ✓；ESP relative RMS = 0.140
- 与文献苯 RESP 电荷（C≈-0.115~-0.15, H≈+0.115~+0.15）一致

### 端到端
- smi2data input.yaml → data.lmp + check_report 通过（电荷/原子数/段计数/类型计数全绿）

## 4. RESP2 实测结果与回退方案

**实测（2026-08-09，AmberTools 26.0 conda-forge）**：
- `antechamber -L` 电荷方法列表：resp / bcc / cm1 / cm2 / esp / mul / gas / abcg2 / rc / wc / dc —— **无 resp2**
- 直接跑 `-c resp2`：`Fatal Error! Unknown charge method (resp2)`
- QUICK 隐式溶剂：`SOLVENT=WATER` 关键词被忽略（带/不带能量完全一致 -37.588557866 a.u.），QUICK 无 PCM/SMD
- RESP2 官方实现路径（Gaussian SCRF(PCM,Water) 输出 → antechamber）需要 Gaussian 许可，本环境没有

**回退方案**：
1. **默认回退 RESP**（气相 HF/6-31G* ESP 两阶段拟合，本报告第 3 节已验证）：电荷质量好（乙醇 RMS 0.203、苯 RMS 0.140），满足阶段 3 验收"仅电荷不同、总和=净电荷"
2. **RESP2 近似**（若博士需要）：装开源 PCM 能力 QM（PSI4/PySCF 的 PCM，或 ORCA 学术免费许可）→ 算水中 HF/6-31G* ESP → smi2data 新增 ESP 源读入接口（.esp 格式已通用，resp 两阶段链路不变）→ 溶剂 RESP 拟合。工作量：一个新 ESP 转换器 + 配置项
3. **RESP2 的非极性 H 特征**：现有 respgen 等效约束已保证 CH3/CH2 的 H 分别等效（乙醇实测），近似满足 RESP2 的合并意图

## 5. 遗留

- RESP2：记录回退方案（见上），待博士决定是否走 PSI4/PySCF/ORCA 路径
- QUICK vdW 表面点数 9777（乙醇）——RESP 拟合足够；大分子性能待观察
