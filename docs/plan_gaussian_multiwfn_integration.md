# 计划：集成 Gaussian 16 + Multiwfn 到 getlmp

> 状态：规划中（2026-08-10）
> 依据：http://sobereva.com/531（RESP2 思想与 Multiwfn 计算）、http://sobereva.com/441（RESP 原理）
> 执行者：换上下文后的 MDToolAgent，按本计划逐步实现

## 1. 背景与目标

**现状**：RESP 电荷用 QUICK（HF/6-31G*，免费）+ antechamber resp 两阶段拟合；
RESP2 不支持（antechamber 无 resp2 方法、QUICK 无 PCM）；ESP 可视化用 esp_cub.py
（点电荷库仑近似）。

**目标**（531 路线）：
- Gaussian 算波函数（B3LYP-D3(BJ)/def2-TZVP，气相 + 可选 PCM 溶剂）
  → Multiwfn 拟合 RESP / RESP2 → 电荷写回 mol2
- Multiwfn 同时承担 ESP 可视化（密度 0.001 闭合等值面 + 表面 ESP → vtx.pdb，
  解决"缺闭合几何来源"的历史问题）
- 保留 QUICK 作为 fallback 引擎（`qm: quick`），不破坏旧路径

## 2. 环境事实（已探测，2026-08-10）

| 项 | 值 |
|---|---|
| 项目目录 | `/home/yryd/packages/getlmp/`（main.py 入口） |
| conda 环境 | `smi2data`（`/home/yryd/packages/miniconda3/envs/smi2data`） |
| CPU | AVX2 支持，8 核，11G 内存，945G 可用磁盘 |
| G16 安装包 | `/mnt/d/BaiduNetdiskDownload/014-gaussian/03-gaussian16linux+09Linux/G16-C01-AVX.tbJ` |
| G16 备选包 | 同目录 `gaosi16-A03-AVX2.tbz`（A03）、`gaosi16-A03-SSE42.tbz`、`g16-b01-sse.tar.gz` |
| G09 包 | 同目录 `gaosi09-C01/D01/E01_EM64T`（不需要，用 G16） |
| Multiwfn | `/mnt/c/Users/youra/Desktop/Multiwfn_2026.7.15_bin_Linux_noGUI.zip`（noGUI 二进制，免编译） |
| Multiwfn 备选 | `/mnt/c/Users/youra/Desktop/Multiwfn_2026.7.15_src_Linux/`（源码，含 Makefile + dislin_d-11.0.a + noGUI/） |

## 3. 阶段 0：前置

1. 确认 git 状态：`cd /home/yryd/packages/getlmp && git status`（历史首提交 c872f75），
   改动前先提交一个 checkpoint。
2. 建软件目录：`mkdir -p /home/yryd/packages/soft/{g16,multiwfn,scratch}`。
3. 确认 conda 可用：`source /home/yryd/packages/miniconda3/etc/profile.d/conda.sh && conda activate smi2data`。

## 4. 阶段 1：安装 Gaussian 16（WSL）

1. 解压：`tar xjf /mnt/d/.../G16-C01-AVX.tbJ -C /home/yryd/packages/soft/`
   （`.tbJ` 即 tbz；若损坏换 A03-AVX2）
2. 环境变量脚本 `scripts/g16env.sh`（供代码 subprocess 复用，不污染全局）：
   ```bash
   export g16root=/home/yryd/packages/soft
   export GAUSS_EXEDIR=$g16root/g16/bsd:$g16root/g16/utility:$g16root/g16
   export GAUSS_SCRDIR=/home/yryd/packages/soft/scratch
   export LD_LIBRARY_PATH=$g16root/g16/bsd:$g16root/g16
   export PATH=$g16root/g16:$PATH
   ```
3. 运行安装脚本：`$g16root/g16/bsd/install`（编译缺失的二进制，需要 gcc）。
4. **license（最大风险）**：确认 `license.dat` 存在且有效。网盘包通常自带；
   若无则 g16 无法运行 → 必须请用户提供合法 license，或回退到 `qm: quick` + Multiwfn。
5. 验证：跑最小 gjf（水分子单点 HF/6-31G*），确认 `Normal termination`；
   验证 `formchk` 可执行（.chk → .fch）。

## 5. 阶段 2：安装 Multiwfn（WSL）

1. 解压：`unzip /mnt/c/Users/youra/Desktop/Multiwfn_2026.7.15_bin_Linux_noGUI.zip -d /home/yryd/packages/soft/multiwfn/`
2. 检查可执行：`chmod +x`、`ldd` 查缺失库（libgomp 等）。
3. 验证：用管道喂输入（如 `7\n18\n1\n0\n`）对测试 fch/molden 跑 RESP，确认输出电荷文件。
4. 若 noGUI 二进制缺库跑不起来 → 备选：用桌面 src_Linux 源码编译（gfortran + noGUI 模式）。
5. 位置约定：config 可配 `multiwfn_path`；默认探测 PATH / `~/packages/soft/multiwfn/`。

## 6. 阶段 3：代码集成（核心）

### 6.1 config.py
新增 `QmCfg` dataclass 与 `Config.qm` 字段，yaml 支持 `qm:` 段：
```python
@dataclass
class QmCfg:
    engine: str = 'gaussian'   # gaussian / quick
    g16root: str = ''          # 空=自动探测（env/默认路径）
    method: str = 'b3lyp'      # 531 推荐 B3LYP-D3(BJ)
    basis: str = 'def2-TZVP'   # 阴离子可配 ma-TZVP
    opt: bool = False          # 是否先几何优化（默认单点）
    solvent: str = ''          # 空=气相；'water'/'ethanol'（RESP2 需要）
    resp2: bool = False        # RESP2（需 solvent 非空）
    delta: float = 0.5         # RESP2 δ
    multiwfn_path: str = ''    # 空=自动探测
```
- `SUPPORTED_CHARGE` 允许 `resp2`（仅当 `qm.engine=gaussian` 且 `qm.resp2=True`、
  `qm.solvent` 非空时通过校验；否则报错说明原因）。
- 校验：`resp` 时 engine 必须 gaussian 或 quick；`bcc/abcg2` 不依赖 qm。

### 6.2 molecule_layer.py：RESP 分支重写
`_build_resp_charges()` 改为按 `qm.engine` 分派：
- **`engine == 'gaussian'`（新主路径）**：
  1. `_write_gaussian_input()`：写 gjf
     （`%nprocshared=8`、`%mem=4GB`、`# B3LYP em=GD3BJ def2-TZVP`
     + `scrf=solvent=xxx`（若溶剂）+ `pop=MK` 或 `iop(6/33=2)` 输出 ESP 点
     + `charge multiplicity`，坐标来自 mol2，不重排原子序）
  2. `_run_gaussian()`：subprocess 带 g16env 环境，检查 `Normal termination`
  3. `_run_formchk()`：.chk → .fch
  4. `_multiwfn_resp()`：Multiwfn 7→18→1（标准两步 RESP），解析电荷
- **`engine == 'quick'`（旧路径保留）**：现有 QUICK 代码原样作为 fallback。
- `_apply_charges()` 不变。

### 6.3 RESP2（新功能）
- 两个单点：气相 fch + 溶剂 fch（PCM）
- Multiwfn 7→18 输入两文件 → RESP2 δ 混合（Multiwfn 支持，参考 calcRESP.sh 逻辑）
- 输出电荷数 = 原子数，守恒校验同 RESP

### 6.4 新模块 multiwfn.py（封装）
- `find_multiwfn()`：探测 config / PATH / 默认路径，友好报错
- `resp_from_fch(fch, out_chg) -> list[float]`
- `resp2_from_fch(fch_gas, fch_solv, delta) -> list[float]`
- `esp_surface_pdb(fch, vtx_pdb)`：密度 0.001 等值面 + 表面 ESP → vtx.pdb
  （B 因子存 ESP，VMD Beta 着色，闭合曲面）
- `esp_cube(fch, cube)`：严格 QM ESP cube（替代点电荷近似，Multiwfn 主功能计算）
- 统一：管道喂 Multiwfn 输入（`7\n18\n...\n0\n`），解析输出，超时保护

### 6.5 ESP 可视化改造（esp 模块）
- 单分子 RESP + `qm=gaussian` 时：Multiwfn 生成 vtx.pdb + cube（严格 ESP）
- 旧 esp_cub.py 点电荷近似保留为 fallback 或删除（实现时确认：若 Multiwfn 路径全覆盖则删）

### 6.6 organize_output / main.py
- 新产物归类：`data.lmp`、mol2、`*.fch`、`*.chg`、`vtx.pdb`、`*.cub` 保留；
  gjf/out/scratch 进 `_others/` 或清理
- main.py 无结构性改动

## 7. 阶段 4：测试与验收

1. **回归**：`pytest tests/`（甲烷/乙醇/混合 3 用例，bcc 路径不受影响，须全过）
2. **RESP 专项**：
   - 乙醇 `charge_method: resp, qm: {engine: gaussian}` → 电荷守恒、
     RMS 输出、与 QUICK resp 结果趋势对比（记录数值）
   - 苯 → 等价性（6C 电荷相同、6H 相同，Multiwfn 自动约束）
3. **RESP2 专项**：乙醇 气相+water → 电荷 = 0.5*(gas+solv)，守恒
4. **ESP 可视化**：vtx.pdb 在 VMD 显示闭合曲面 + Beta 着色（人工确认）
5. **验收标准**：
   - data.lmp 电荷和 = 体系净电荷（容差 1e-4）；原子数守恒；LAMMPS 可读
   - RESP/RESP2 电荷物理合理（对照 531 示例与 QUICK 结果）
   - `qm: quick` 旧路径零回归

## 8. 阶段 5：文档与归档

- `docs/dev_notes.md`：更新 RESP 章节（Gaussian+Multiwfn 流程、531 依据、
  RESP2 说明、环境变量脚本用法）
- `README.md`：charge_method 支持矩阵（resp2 可选）、`qm:` 配置示例
- `examples/`：新增 `gaussian_resp.yaml`、`resp2.yaml`
- 检查报告归档到项目目录，git 提交

## 9. 风险与回退

| 风险 | 影响 | 回退 |
|---|---|---|
| **Gaussian license 缺失/无效**（最大） | g16 无法运行 | `qm: quick` + Multiwfn（Multiwfn 读 molden 拟合 RESP）仍可交付；或请用户提供 license |
| G16 C01 AVX 包损坏/WSL 不兼容 | 装不上 | 换 A03-AVX2 / SSE42 版 |
| Multiwfn noGUI 缺动态库 | 跑不起来 | 桌面源码编译（gfortran + noGUI）；或改用 bin_Win64 + wine（不推荐） |
| 环境变量污染 | 影响其他软件 | 统一走 `scripts/g16env.sh`，仅 subprocess 内注入 |
| 内存不足（大分子） | g16 OOM | 小分子（<50 原子）11G 内存足够；大分子提示降低核数/基组 |

## 10. 交付物清单

- [ ] g16 + Multiwfn 在 WSL 可运行（含验证记录）
- [ ] config 支持 `qm:`（engine/g16root/method/basis/solvent/resp2/delta/multiwfn_path）
- [ ] molecule_layer RESP 分支支持 Gaussian + Multiwfn 流程（QUICK 保留）
- [ ] multiwfn.py 封装（RESP/RESP2/ESP 曲面/ESP cube）
- [ ] ESP 可视化产出 vtx.pdb（闭合曲面）
- [ ] 测试：回归 3/3 + RESP/RESP2 专项通过
- [ ] 检查报告 + dev_notes/README/examples 更新 + git 提交
