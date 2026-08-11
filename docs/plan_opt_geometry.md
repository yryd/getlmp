# 计划：qm.opt 优化几何回写 mol2（最终 data.lmp 用 QM 优化结构）

> 状态：规划中（2026-08-11）
> 需求：`qm.opt: true` 时不仅用优化几何算 RESP 电荷，还要把优化后的坐标
> **回写 mol2**，使最终 `data.lmp` 的原子坐标 = QM 优化结果（当前只优化电荷拟合，
> 坐标仍是 RDKit 构象，opt 计算浪费）。
> 默认仍为 `false`（单点），开关语义不变。

## 1. 现状（代码事实）

| 位置 | 现状 |
|---|---|
| `src/config.py:63` | `QmCfg.opt: bool = False` 字段（默认关）✅ |
| `src/config.py:124` | yaml 解析 `cfg.qm.opt` ✅ |
| `src/molecule_layer.py:323-324` | `_write_gaussian_input`：`if qm.opt: route += ' opt'` ✅ |
| `src/molecule_layer.py:111` | `_apply_charges(mol2, charges)`：**只写电荷、坐标不变** ← 回写接入点 |
| `src/molecule_layer.py:248` | `_mol2_atoms(mol2)`：已有 mol2 ATOM 段读取 |
| `_build_resp_charges_gaussian` (424) | `gjf → g16 → formchk → RESP`；插入点：`_run_gaussian` 后 |
| `_build_resp2_charges` (452) | gas + solv 两个 gjf 各自跑；opt 时两者都优化 |
| QUICK 路径 | 无几何优化器，opt 不支持（保持现状） |

**当前行为**：`opt: true` 时 g16 优化几何 → 更准的 ESP → 更准的电荷；
但 mol2/data.lmp 坐标仍是 RDKit ETKDG 构象（`_apply_charges` 不碰坐标）。

## 2. 方案设计

### 2.1 坐标来源：g16 log 的 "Input orientation" 最后一块
- g16 优化每步输出一个 `Input orientation`（单位 Å，含原子序 + 坐标）；
  取**最后一块**即优化收敛几何。
- **必须加 `nosymm`**：g16 默认对称性会把原子重排（苯 6H 会换序），
  导致坐标与 mol2 原子顺序错位。`nosymm` 后 Input orientation 与输入顺序一致。
  代价：对称分子优化略慢（小分子可忽略）。
- 备选（不采用）：`.fch` 的 `Current cartesian coordinates`（Bohr、按原子序排，
  仍有顺序问题且多一次单位转换）。

### 2.2 优化完成判断
- 现有 `_run_gaussian` 检查 `Normal termination`——但**优化不收敛时 g16 也可能
  正常退出**（报 `Optimization stopped` / 无 `Stationary point found`）。
- 新增：解析 log 出现 `Stationary point found` 才算优化成功；否则报错
  （给出 log 尾部，提示调小 `opt` 步长/检查初始结构）。

### 2.3 写回 mol2
- 新函数 `_update_mol2_coords(mol2, xyz)`：与 `_apply_charges` 对称——
  解析 ATOM 段替换第 3–5 列坐标（保留类型/电荷/键），断言原子数一致。
- 调用时机：**先回写坐标，再做 formchk/RESP**（电荷拟合用优化后的 fch，
  坐标也用优化后的——两者一致，避免"优化电荷配旧坐标"）。

### 2.4 RESP2 双单点：回写哪个坐标？
- gas 与 solv 各自优化到不同几何。**默认回写溶剂（PCM）优化坐标**
  （更贴近真实环境，RESP2 官方惯例）；gas 单点仅用于电荷混合。
- 实现：RESP2 流程按 gas→solv 顺序跑，取**最后一个** gjf（=溶剂）的优化坐标回写。
- 不引入新配置（MVP 不过度设计）。

### 2.5 QUICK 路径
- 无优化器，`opt` 保持忽略；加一行提示（`[info] quick 引擎不支持 opt，已忽略`）。

## 3. 代码改动清单（估算 ~90 行新增）

| 文件 | 改动 | 行数 |
|---|---|---|
| `src/molecule_layer.py` | 新 `_parse_g16_opt_coords(log)`：扫描 Input orientation 块，取最后一块坐标 | ~40 |
| `src/molecule_layer.py` | 新 `_update_mol2_coords(mol2, xyz)`：ATOM 段坐标写回 | ~30 |
| `src/molecule_layer.py` | `_write_gaussian_input`：`if qm.opt: route += ' opt nosymm'` | ~1 |
| `src/molecule_layer.py` | `_build_resp_charges_gaussian`：`_run_gaussian` 后插入「解析→回写」+ 收敛检查 | ~10 |
| `src/molecule_layer.py` | `_build_resp2_charges`：solv 单点后回写（取最后一个 gjf 坐标） | ~8 |
| `src/molecule_layer.py` | `_build_resp_charges_quick`：opt 忽略提示 | ~2 |
| `docs/yaml_config.md` | `opt` 键说明补「回写 mol2」+ RESP2 回写溶剂 | ~5 |
| `docs/dev_notes.md` | 4.1 补 opt 回写说明与验证数据 | ~10 |
| `docs/install.md` | 无改动 | 0 |

## 4. 验证方案

1. **opt:false 零回归**：`tests/test_build_systems.py` 3/3 + 乙醇 RESP（坐标应与
   之前一致，电荷数值一致）
2. **乙醇 opt:true**：跑通全链路；断言 mol2/data.lmp 坐标 ≠ RDKit 原始坐标
   （优化后 C-C/C-O 键长更接近实验值）；电荷守恒；tleap 正常
3. **苯 opt:true（对称分子坑）**：nosymm 后原子顺序不乱；写回坐标后 RESP 电荷
   仍 6C/6H 等价（≈±0.109）；tleap 无原子错位报错
4. **RESP2 opt:true**：双单点都优化；mol2 坐标 = 溶剂优化坐标；δ 混合电荷守恒
5. **边界**：构造一个不收敛场景（可选，至少确认报错信息存在）
6. 检查报告：`work/_others/check_report.txt` 记录 opt 状态与坐标来源

## 5. 工作量评估

| 项 | 估算 |
|---|---|
| 代码实现（molecule_layer.py） | ~90 行，0.5–1 小时 |
| 乙醇 + 苯 + RESP2 opt 冒烟 | 每次 ~2–3 分钟（B3LYP/def2TZVP 小分子），含排错 ~1 小时 |
| 文档 + 示例 + 提交 | ~0.5 小时 |
| **合计** | **约 2–3 小时**（含验证与回归） |

难度：**低**。风险集中在 g16 log 解析健壮性与对称分子原子顺序（nosymm 已规避）。

## 6. 风险与回退

| 风险 | 影响 | 处理 |
|---|---|---|
| g16 log 格式差异（多版本） | 坐标解析失败 | 本机 G16 RevC.01 实测；解析失败报错含 log 尾部，不静默 |
| 对称分子原子重排 | 坐标错位 | `nosymm` 规避；苯专项验证 |
| opt 不收敛 | 无 Stationary point | 报错提示（现有 Normal termination 检查之外补强） |
| nosymm 拖慢优化 | 时间增加 | 小分子可忽略；文档注明大分子可权衡 |
| 回写坐标破坏 tleap | 拓扑失败 | 原子数断言 + tleap 报错天然兜底；苯专项验证 |

## 7. 交付物

- 代码：molecule_layer.py（坐标解析/写回/编排）
- 文档：yaml_config.md（opt 说明）、dev_notes.md（4.1 验证数据）
- 验证：乙醇/苯/RESP2 opt 冒烟 + 回归 3/3 + 检查报告
- 提交信息含工作量与验证记录
