# getlmp YAML 配置参考

> 唯一入口：`python main.py input.yaml`。本文列出**全部可用配置项**：类型、默认值、
> 可选值、作用与组合约束。来源：`src/config.py`（`load_config` 校验逻辑）。
> 更新日期：2026-08-14。

---

## 0. 速查：最小示例与完整示例

**最小示例**（单分子，AM1-BCC）：

```yaml
name: ethanol
forcefield: gaff2
charge_method: bcc
molecules:
  - smiles: CCO
    name: ethanol
```

**溶剂体系示例**（溶质 + TIP3P 水 + NaCl）：

```yaml
name: pip_water
forcefield: gaff2
charge_method: bcc
molecules:
  - smiles: C1CNCCN1
    name: PIP
    count: 20
water:
  model: tip3p       # tip3p / spce / opc3（二期: opc / tip4pew / tip4pd）
  count: 800
ions:
  - name: Na+
    count: 20
  - name: Cl-
    count: 20
packmol:
  density: 1.0       # 水+溶质密度 1.0 g/cm³ → 自动算盒
workdir: work
output: data.lmp
```

**完整示例**（RESP2，展示全部常用段）：

```yaml
name: resp2_ethanol
forcefield: gaff2
charge_method: resp2
net_charge: 0
qm:
  engine: gaussian
  method: b3lyp
  basis: def2TZVP
  solvent: water
  resp2: true
  delta: 0.5
esp:
  enabled: true
  spacing: 0.3
  buffer: 1.5
molecules:
  - smiles: CCO
    name: ethanol
    count: 1
output: data.lmp
workdir: work
seed: 2026
organize_output: true
```

---

## 1. 顶层选项

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `name` | str | `system` | 体系名。用于主产物命名（`{name}.fch`/`{name}.chg`、检查报告标题等） |
| `forcefield` | str | `gaff2` | 力场：`gaff2` / `gaff` / `reaxff` |
| `charge_method` | str | `bcc` | 电荷方法：`bcc`(AM1-BCC) / `abcg2` / `resp` / `resp2` / `none`（`reaxff` 时必须为 `none`） |
| `net_charge` | int | `0` | 分子净电荷（传给 antechamber `-nc` 与 g16 输入；多分子体系按每个分子拷贝计） |
| `molecules` | list | **必填**（除非配了 `water`/`ions`） | 分子列表，见 §2（纯溶剂体系可省略） |
| `output` | str | `data.lmp` | 输出的 LAMMPS data 文件名（相对 `workdir`） |
| `workdir` | str | `work` | 工作目录（相对配置文件所在目录） |
| `seed` | int | `2026` | 随机种子：RDKit ETKDG 3D 构象 + packmol 装盒共用 |
| `organize_output` | bool | `true` | 跑完自动整理：根目录只留主产物（data.lmp、mol2、system.xyz、ESP 正式产物、fch/chg），其余全部移入 `workdir/_others/` |
| `organize_backup` | bool | `false` | 与 `_others/` 已有文件同名时：`false`=覆盖只留最新；`true`=追加 HHMMSS 时间戳保留历史多份 |
| `reuse_molecule` | bool | `false` | 分子层复用：`workdir` 已有**同配置指纹**的 mol2/frcmod/波函数时跳过重算（antechamber/QM 不重跑）。改 SMILES/力场/电荷方法/净电荷/seed/QM 设置后指纹不匹配会自动重算。默认关 |
| `buffer` | float | `3.8` | 单分子盒 padding（Å）。仅 `count=1` 时使用（`[坐标极值 ± buffer]` 推算盒）；多分子由 `packmol.box` 决定 |
| `reax_elements` | list[str] | `[C, H, O, N, S, P, F, Cl, Br, I]` | ReaxFF 专用：data 类型号→元素顺序，**必须与你的 `ffield.reax` 元素顺序一致**；可覆盖 |
| `reax_atom_style` | str | `charge` | ReaxFF 专用：`charge`=Atoms 段 6 列（无 mol-id，需分子分组时配 `full`）；`full`=7 列带 mol-id |
| `qm` | dict | `{}` | 量子化学引擎段（RESP/RESP2），见 §3 |
| `esp` | dict | `{}` | ESP 可视化导出段，见 §4 |
| `packmol` | dict | `{}` | 多分子装盒段，见 §5 |

### 校验规则（顶层）
- `forcefield` / `charge_method` 取值非法直接报错；
- `forcefield: reaxff` 时 `charge_method` 必须 `none`（QEq 在模拟中计算）；
- 缺少 `molecules` 段时，必须有 `water` 或 `ions`（纯溶剂体系），否则报错。

---

## 2. molecules 段（分子列表，必填）

每个条目一个字典：

```yaml
molecules:
  - smiles: CCO        # SMILES（RDKit 解析）
    name: ethanol      # 分子名（默认 mol1/mol2...）
    count: 1           # 拷贝数（>1 自动启用 packmol 装盒）
    resname: ETH       # 残基名（antechamber -rn；默认取 name 前 3 字符大写）
  - xyz: mol.xyz       # 或直接给 xyz 坐标文件（与 smiles 二选一）
    name: PIP
    count: 5
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `smiles` | str | 与 `xyz` 二选一 | SMILES 字符串，RDKit 解析并生成 3D 构象 |
| `xyz` | str | 与 `smiles` 二选一 | **xyz 坐标文件路径**（相对 yaml 所在目录或绝对路径）。元素符号 + 坐标即可（第 2 行注释可选），RDKit 读入后走与 SMILES 完全相同的分子层链路（加氢不适用——xyz 原子数即最终原子数；antechamber 从坐标建拓扑）。**适用于：已有构象/来自其他软件的体系、或 SMILES 无法表达的坐标** |
| `name` | str | `mol{序号}` | 分子名（也用于中间文件命名：`{name}.sdf/.mol2/.fch` 等） |
| `count` | int | `1` | 分子拷贝数。**任一条目 count>1 即视为多分子体系**（自动 packmol，需配 `packmol.box`） |
| `resname` | str | `name[:3].upper()` | 残基名，antechamber `-rn` 用；默认取 name 前 3 字符大写 |

> `smiles` 与 `xyz` 二选一，同时给会报错。xyz 输入不支持无氢（如金属原子）自动加氢——xyz 里的原子就是最终原子。

---

## 2.5 water 段（水溶剂模板，可选）

> 出现即视为**多分子体系**（与 molecules 一起 packmol 装盒）。也可只有 `water` 无
> `molecules`（纯水体系）。水参数不自己做——**直接用 AmberTools 内置水模型**
> （`leaprc.water.*` + box off 库），与 tleap 加载完全一致，保证 prmtop→data.lmp 自洽。

```yaml
water:
  model: tip3p    # 一期可用: tip3p / spce / opc3；二期预留: opc / tip4pew / tip4pd
  count: 3000     # 水分子数
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `model` | str | **必填** | 水模型。**一期（3-site）**：`tip3p`（经典默认）/ `spce` / `opc3`；**二期预留（4-site）**：`opc` / `tip4pew` / `tip4pd`——填了会直接报错提示（4-site 含虚拟位点 EP，LAMMPS data 导出需专门链路，二期实现） |
| `count` | int | **必填** | 水分子数 |

**实现要点**（对模拟的影响）：
- 水模型由 tleap 的 `leaprc.water.{model}` 加载（键长/角/电荷/LJ 全来自 AmberTools）；
- 水坐标从对应 box off 库解析（packmol 装盒用），与 tleap 加载的库严格一致；
- 水分子无 O-H 键/角定义（刚性模型 + SETTLE），导出 data.lmp 时**自动补**键角：
  TIP3P 角 k=100 θ=104.52°；SPC/E θ=109.47°；OPC3 θ=109.43°（各模型文献标准值，
  跑刚性水模拟用 `fix rigid/settle` 时角参数被忽略）；
- 数据流：`packmol 装水（模板坐标）→ loadpdb → tleap 按水模板分配类型/电荷 → data.lmp`，
  电荷自动平衡（每分子 O=-2qH），无需 `net_charge` 参与。

---

## 2.6 ions 段（离子溶剂，可选）

```yaml
ions:
  - name: Na+      # atomic_ions.lib 类型名（单原子离子）
    count: 20
  - name: Cl-
    count: 20
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `name` | str | **必填** | 离子名，**必须匹配 AmberTools `atomic_ions.lib` 的类型名**（单原子离子）：`Na+` `K+` `Rb+` `Cs+` `Li+` `F-` `Cl-` `Br-` `I-` `Mg2+` `Ca2+` `Zn2+` `Cu2+` `Fe2+` `Fe3+` `Al3+` 等（运行时动态扫描该库，以库内为准）。**多原子离子**（H3O+、NH4+、SO4²⁻ 等）**不做特殊处理**——用 SMILES 走普通分子层（硫酸根已验证：SMILES 输入 → GAFF2 参数 + 电荷守恒） |
| `count` | int | **必填** | 离子数 |

**参数来源**（重要）：
- 离子坐标：`atomic_ions.lib` 单原子模板（packmol 装盒）；
- **LJ 参数随 tleap leaprc 自动加载**：配 TIP3P 水时是 **Joung-Cheatham 12-6**
  （`frcmod.ionsjc_tip3p`，Na+ ε=0.0874 kcal/mol σ=2.439 Å，Cl- ε=0.0356 σ=4.478 Å，
  amber 默认）；配 SPC/E/OPC3 水时是各自配套的离子参数（Li/Merz 12-6 等）。
  检查报告会列出 data.lmp 中离子的实际 LJ 参数供核对。
- **12-6-4 说明**：Amber 离子参数族里有 12-6-4 变体（含 C4 项），但 LAMMPS 无原生
  12-6-4 支持（需 NBFIX 类补丁，有坑），本项目默认走标准 12-6 JC/Li-Merz 链路。

---

## 3. qm 段（量子化学引擎，RESP/RESP2 用）

> 仅 `charge_method: resp` / `resp2` 时生效；`bcc` / `abcg2` 不依赖 qm（bcc 用 sqm 半经验，
> abcg2 用内置查表）。**不写 qm 段 = 全部默认**（engine=gaussian）。

```yaml
qm:
  engine: gaussian     # gaussian（默认，G16+Multiwfn）/ quick（QUICK 旧路径回退）
  g16root: ''          # G16 根目录（空=自动探测；通常不用配）
  method: b3lyp        # QM 方法（531 推荐 B3LYP-D3(BJ)）
  basis: def2TZVP      # 基组（G16 命名；文献写法 def2-TZVP 会自动转）
  opt: false           # true=先几何优化再算 ESP，且优化坐标回写 mol2（最终 data.lmp 用 QM 优化结构）；默认 false=单点
  solvent: ''          # 空=气相；'water'/'ethanol'（PCM 隐式溶剂；RESP2 必需）
  resp2: false         # RESP2 开关（需 solvent 非空）
  delta: 0.5           # RESP2 δ 混合系数 [0,1]
  multiwfn_path: ''    # Multiwfn_noGUI 路径（空=自动探测）
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `engine` | str | `gaussian` | `gaussian`：G16 单点（`pop=MK IOp(6/33=2,6/42=6)`）+ Multiwfn 拟合，531 路线主路径；`quick`：QUICK HF/6-31G* + resp 两阶段，旧路径回退（免费、不依赖 G16，但**无 RESP2、无隐式溶剂**） |
| `g16root` | str | `''` | G16 安装根目录。**一般留空**：代码按 PATH 找 `g16`（见 `docs/install.md`），只需 G16 已在 PATH 中。仅非常规安装位置（不在 PATH）时才显式指定 |
| `method` | str | `b3lyp` | 泛函/方法名，G16 route 直接使用。531 路线推荐 `b3lyp`（代码自动补 `em=GD3BJ` 色散校正） |
| `basis` | str | `def2TZVP` | 基组。注意 **G16 不认带连字符的 `def2-TZVP`**，代码会自动把 `def2-` 转 `def2`（兼容文献写法）。阴离子体系可配 `ma-def2TZVP` 等 |
| `opt` | bool | `false` | `true` 时 g16 先做几何优化（route 加 `opt nosymm`），**优化坐标回写 mol2**，最终 `data.lmp` 原子坐标 = QM 优化结果（不再浪费优化）。默认 `false`=单点（几何=RDKit ETKDG 构象→mol2，不经 QM 优化）。RESP 惯例为单点，一般保持默认；需要"QM 优化结构做最终建模"时开启。注意：`nosymm` 保证原子顺序与 mol2 一致（苯等对称分子不会重排）；RESP2 开启时回写**溶剂(PCM)优化坐标**（更贴近真实环境）；`engine: quick` 无优化器，该键被忽略 |
| `solvent` | str | `''` | 隐式溶剂名（G16 SMD 关键词，如 `water`/`ethanol`）。非空时 route 加 `scrf=(smd,solvent=xxx)`。**RESP2 必需**（溶剂单点） |
| `resp2` | bool | `false` | RESP2 开关（与 `charge_method: resp2` 配套，见 §6 校验） |
| `delta` | float | `0.5` | RESP2 混合系数：`q = (1-δ)·q_gas + δ·q_solv`。取值范围 [0,1]，越接近 1 越偏向溶剂化 |
| `multiwfn_path` | str | `''` | Multiwfn_noGUI 可执行文件路径。**一般留空**：代码按 PATH 找（见 `docs/install.md`）。仅不在 PATH 时才显式指定 |

---

## 4. esp 段（ESP 可视化导出）

> 仅**单分子**（全部 count=1）且 `charge_method: resp/resp2` 时生效；多分子自动跳过（打印提示）。

```yaml
esp:
  enabled: true    # 不需要可视化时 false（不影响 RESP 电荷拟合本身）
  pt: false        # pt 法产物（mol.pdb/vtx.pdb），默认关（基本不用；需要时开）
  spacing: auto    # auto=按原子数分档（≤20:0.25 / 21-40:0.3 / >40:0.4）| 0.15~0.8 显式
  timeout: auto    # auto=按原子数分档（≤20:600 / 21-40:1800 / >40:3600）| 秒数；0=不限
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `enabled` | bool | `true` | 是否导出 ESP 可视化产物。`false` 只关可视化，**不影响电荷拟合** |
| `pt` | bool | `false` | 是否导出 pt 法产物（默认关）。pt 法输出 `pt/mol.pdb` + `pt/vtx.pdb` |
| `spacing` | `auto`/float | `auto` | iso/pt 法**表面格点间距 Å**。`auto`：按原子数分档（≤20→0.25，21-40→0.3，>40→0.4）；显式：0.15~0.8（Multiwfn 手册建议 0.15~0.25 精细，0.4~0.6 大分子加速）。2026-08-14：Multiwfn 默认 0.25 对 ≥1000 基函数大分子表面点数爆炸会超时，故引入分档 |
| `timeout` | `auto`/int | `auto` | Multiwfn 单次调用超时**秒**。`auto`：按原子数分档（≤20→600，21-40→1800，>40→3600）；`0`=不限（不推荐，卡死无兜底） |

**产物**（统一 Multiwfn 导出，engine=gaussian 用 `.fch`、quick 用 `.molden`），位于 `_others/electrostatic_potential/`：
- `iso/density.cub` + `iso/esp.cub`：iso 法——严格 QM 电子密度与 ESP 网格（同一网格，VMD density 等值面 0.001 + esp 着色）。**默认生成**。
- `pt/mol.pdb` + `pt/vtx.pdb`：pt 法——分子结构与闭合等值面顶点（B 因子=ESP kcal/mol，VMD Beta 着色）。**需 `esp.pt: true`**。

> 旧字段 `buffer` 已废弃移除；yaml 里遗留的 `buffer` 键会被忽略。

---

## 5. packmol 段（多分子装盒）

> `count>1` 自动启用（`enabled` 无需手写）；启用时**必须**给 `box` 或 `density`（二选一，`density` 优先）。

```yaml
packmol:
  enabled: true     # 多分子自动开；单分子想强制装盒可手写 true
  preset: bulk      # 当前仅支持 bulk（slab/interface 后续阶段）
  box: [0, 0, 0, 60, 60, 60]   # [xlo, ylo, zlo, xhi, yhi, zhi] Å（与 density 二选一）
  density: 1.0      # 目标密度 g/cm³（>0 时自动按 总质量/密度 算立方盒，忽略 box）
  seed: 2026        # 装盒随机种子（默认取顶层 seed）
  tolerance: 2.0    # 装盒容差（Å，分子间最小间距）
  nloop0: 1000      # packmol 初始随机放置尝试次数（默认 20 对高密度大体系常不够）
  inp_file: work/my.inp   # 自定义 packmol inp（可选；非空时跳过自动生成，直接用该文件）
```

| 键 | 类型 | 默认 | 说明 |
|----|------|------|------|
| `enabled` | bool | `false`（任一条目 count>1 时自动 true） | 是否 packmol 装盒。单分子体系也可显式 `true` 强制装盒 |
| `preset` | str | `bulk` | 装盒模式。**当前仅 `bulk`**；`slab`/`interface` 在规划中（填了会报错） |
| `box` | list[6] | 无 | 盒尺寸 `[xlo, ylo, zlo, xhi, yhi, zhi]`（Å）。与 `density` 二选一（`density` 优先） |
| `density` | float | `0.0` | 目标密度 g/cm³。>0 时按 **L = (总质量/密度)^(1/3)** 自动算立方盒并写回 `box`（write_inp 与 data.lmp 盒边界共用），质量按各分子 mol2 元素组成 × count 求和 / N_A（ReaxFF 按 `reax_elements` 元素表）。装不下（packmol 报错）时调小 `density`（盒子变大）或调小 `tolerance` |
| `seed` | int | 顶层 `seed` | 装盒随机种子 |
| `tolerance` | float | `2.0` | packmol tolerance（Å，分子间最小间距） |
| `nloop0` | int | `1000` | packmol 初始随机放置循环数（官方默认 20，高密度/大分子体系常不足，报错会提示加大） |
| `inp_file` | str | `''` | 自定义 packmol inp 路径（相对配置文件所在目录或绝对路径）。非空时**跳过自动生成** write_inp，直接用该文件跑 packmol，适合加约束（如 `fixed` 固定分子位置）、改 number、同一类型拆多个 structure 块。**规则**：structure 文件名为 `{name}.xyz`（与 molecules 的 name 对应）；number 以 inp 内为准（可不同于 yaml count）；盒边界以 yaml `box` 为准（inp 里改 inside box 不会同步到 data.lmp）；**与 `density` 互斥** |

---

## 6. 组合约束（load_config 校验）

| 组合 | 约束 |
|------|------|
| `charge_method: resp2` | 必须同时 `qm.engine: gaussian`（QUICK 无 resp2）、`qm.resp2: true`、`qm.solvent` 非空 |
| `forcefield: reaxff` | `charge_method` 必须 `none`；ReaxFF 路径忽略 qm/esp 段 |
| `packmol.enabled`（多分子） | 必须 `packmol.box` 或 `packmol.density`（>0）；`preset` 只能 `bulk` |
| `packmol.inp_file` | 仅多分子可用（单分子报错）；文件必须存在；structure 文件名须为 `{name}.xyz`（对应 molecules 的 name）；**与 `packmol.density` 互斥** |
| `qm.delta` | 必须在 [0,1] |
| `esp.spacing` | `auto` 或 0.15~0.8（Å） |
| `esp.timeout` | `auto` 或 ≥0 整数（0=不限） |
| `reax_elements` | 字符串列表；重复元素自动去重 |
| `reax_atom_style` | 仅 `charge` / `full` |
| `packmol.box` | 必须 6 个数 |
| `molecules[].smiles/xyz` | 二选一必填（同时给或都不给报错） |
| `molecules` 段 | 可省略——当且仅当配置了 `water` 或 `ions`（纯溶剂/离子体系） |
| `water.model` | 一期 `tip3p`/`spce`/`opc3`；二期预留 `opc`/`tip4pew`/`tip4pd`（直接报错提示） |
| `ions[].name` | 必须是 `atomic_ions.lib` 中的单原子类型（运行时扫描校验） |
| `water` / `ions` 出现 | 强制多分子路径（packmol 装盒 + tleap loadpdb），即使所有 count=1 |

## 7. 力场 × 电荷可用性矩阵

| forcefield | charge_method | 引擎 | 状态 |
|---|---|---|---|
| gaff2 / gaff | `bcc` | sqm（AM1-BCC，半经验） | ✅ 默认 |
| gaff2 / gaff | `abcg2` | antechamber 内置查表 | ✅ |
| gaff2 / gaff | `resp` | G16+Multiwfn（默认）或 QUICK 回退 | ✅ |
| gaff2 / gaff | `resp2` | G16 气相+溶剂双单点 + Multiwfn | ✅ |
| reaxff | `none` | 无（QEq 模拟中算） | ✅ |

---

## 8. 常见注意

1. **相对路径**以配置文件所在目录为基准（`workdir`、`output`）。
2. **`qm.opt: true` 会显著变慢**（几何优化通常需多轮 SCF）；RESP 常规用单点即可。开启后优化坐标会回写 mol2（最终 `data.lmp` 用 QM 优化结构），苯等对称分子已验证原子顺序正确（`nosymm` 规避重排）。
3. **`organize_output: false`** 可关闭自动整理，全部中间文件留在 `workdir/` 便于排查。
4. **`organize_backup: true`** 会在 `_others/` 同名文件追加时间戳保留多份（默认覆盖只留最新）。
5. 多分子体系（任一条目 count>1）**不会导出 ESP 可视化**（当前仅单分子支持）。
6. G16 未安装/未配置环境时，RESP 请用 `qm.engine: quick`；RESP2 必须 gaussian。
7. **自定义 packmol inp 流程**：先不带 `inp_file` 跑一次（`workdir/packmol.inp` 会保留在根目录作模板）→ 复制改名（如 `work/my.inp`）自由修改（加 `fixed` 固定某分子、改 number、同类型拆多块）→ yaml 加 `inp_file: work/my.inp` 再跑。改 inp 后 `packmol.inp`/`my.inp` 都不会被 organize 移走。
8. **`reuse_molecule: true` 复用范围**：仅分子层（mol2/frcmod/波函数）按指纹跳过；**体系层 packmol/tleap 每次都重跑**（改 inp、box、count 不影响分子层复用）。指纹文件 `{name}.fingerprint.json` 与 frcmod 会保留在 workdir 根目录（organize 不移动）。
9. **水/离子体系是"自动多分子"**：yaml 出现 `water`/`ions` 段即走 packmol + tleap loadpdb 链路（即使所有 count=1），必须配 `packmol.box` 或 `packmol.density`。**密度是按"全部分子（溶质+水+离子）总质量/体积"算的**，所以 `density: 1.0` 对"溶质+水"体系即目标密度 1 g/cm³。
10. **水/离子体系不要给 `net_charge`**：水与离子电荷由 tleap 模板自动平衡（水 O=-2qH、离子 ±1/±2），导出校验只检查总电荷≈0。溶质自身若带电（如硫酸根 -2），靠 `ions` 里配反离子中和（如 2×Na+）。
