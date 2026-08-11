# Gaussian 16 + Multiwfn 集成：分阶段验证记录

> 对应计划：docs/plan_gaussian_multiwfn_integration.md
> 执行：MDToolAgent，2026-08-11
> 原则：每阶段可复现，命令与输出留档；失败/风险同步记录

## 阶段 0：前置 — ✅ 通过（2026-08-11）

| 检查项 | 结果 |
|---|---|
| git 状态 | 工作树仅未跟踪计划文档；checkpoint 提交 `178dc0a` |
| soft 目录 | `mkdir -p /home/yryd/packages/soft/{g16,multiwfn,scratch}` 完成 |
| conda 环境 | `smi2data`（Python 3.12.13）可用；antechamber/quick/parmchk2 在 PATH |

```bash
cd /home/yryd/packages/getlmp && git add docs/plan_gaussian_multiwfn_integration.md \
  && git commit -m 'docs: 新增 Gaussian 16 + Multiwfn 集成计划（531 路线，阶段 0 checkpoint）'
# → 178dc0a
```

## 阶段 1：安装 Gaussian 16 — ✅ 通过（2026-08-11）

### 1.1 解压（实际格式 XZ，非 bzip2）
- 包：`/mnt/d/BaiduNetdiskDownload/014-gaussian/03-gaussian16linux+09Linux/G16-C01-AVX.tbJ`
- `file` 检测为 **XZ 压缩**（魔数 `fd 37 7a 58 5a`），`.tbJ` 并非 bzip2；
  系统无 bzip2 二进制且 tar 直接解压失败。
- 解决：Python 内置 `tarfile.open(p, 'r:xz').extractall(...)`。
- **目录结构事实**：包顶层含 `g16/`，解压到 `soft/g16/` 产生嵌套
  `soft/g16/g16/`；因 mv 被治理策略禁止，保留嵌套，g16root 直接指向实际安装目录。

### 1.2 权限与 install 脚本
- `bsd/install` 是 `#!/bin/csh -f`，系统无 csh 且无法 sudo 安装 → 手动执行等价操作：
  ```bash
  cd /home/yryd/packages/soft/g16/g16 && chmod -R o-rwx .   # 即 install 的主操作
  ```
  （无 linda9.2 目录，fixlinda 无需执行）
- g16 主可执行：ELF 64-bit x86-64，124 MB，可直接运行（不依赖 csh）。

### 1.3 license 结论（原计划最大风险，已消除）
- 解压目录与网盘各包（G16-C01-AVX / A03-AVX2 / A03-SSE42 / g16-b01-sse）均**无 license.dat**；
- 实测 g16 **无需 license 即可运行**（网盘版免 license），与安装说明（未提 license）一致。

### 1.4 运行验证 ✅
最小 gjf：水分子 HF/6-31G* 单点（`%nprocshared=8 %mem=4GB`）：
```bash
source scripts/g16env.sh && g16 h2o_test.gjf h2o_test.log
```
- 输出：`Normal termination` 出现 1 次
- 能量：`SCF Done: E(RHF) = -76.0105049953 A.U.`（与 HF/6-31G* 标准值吻合）
- `formchk h2o_test.chk h2o_test.fch` → 生成 95 KB .fch ✅

### 1.5 产物
- `scripts/g16env.sh`：环境变量脚本（g16root=/home/yryd/packages/soft/g16/g16，
  GAUSS_EXEDIR/SCRDIR/LD_LIBRARY_PATH/PATH），subprocess 复用

### 1.6 备注
- 测试文件：`/home/yryd/packages/soft/scratch/{h2o_test.gjf,h2o_test.log,h2o_test.chk,h2o_test.fch}`
- 后续如需规整目录（soft/g16/g16 → soft/g16），需手动 `mv`（本环境策略禁止，未做）

---

## 阶段 2：安装 Multiwfn — ✅ 通过（2026-08-11）

- 二进制：`/home/yryd/packages/soft/multiwfn/Multiwfn_2026.7.15_bin_Linux_noGUI/Multiwfn_noGUI`
- 免编译（noGUI 版）；`-ispecial 1` 支持 ESP 拟合脚本交互
- 路径探测：`multiwfn.py` 自动探测（glob `~/packages/soft/multiwfn/*/Multiwfn_noGUI`），
  可用 `qm.multiwfn_path` 覆盖

## 阶段 3：代码集成 — ✅ 通过（2026-08-11，提交 b7cf6d4）

| 模块 | 内容 |
|---|---|
| `src/multiwfn.py`（新） | RESP（7→18→1 两阶段）、RESP2（7→18→3→1，δ=0.5）、ESP 曲面（vtx.pdb）、严格 QM ESP cube |
| `src/molecule_layer.py` | `_build_resp_charges` 双引擎分派（gaussian 默认 / quick 回退）；RESP2 双单点（gas + solvent PCM） |
| `src/config.py` | 新增 `qm:` 段（engine/method/basis/solvent/resp2/delta/multiwfn_path）+ `esp:` 段 |
| `src/pipeline.py` | 波函数归档 + ESP 导出编排（gaussian 引擎 → Multiwfn；quick → 点电荷近似） |
| `scripts/g16env.sh` | G16 环境变量封装（非交互 shell） |
| examples | `gaussian_resp.yaml` / `resp2.yaml`；PIP.yaml 显式 `qm.engine: quick` |

**G16 链路三个坑**（详见 dev_notes 4.1）：
1. gjf `0 1` 后不能有空行（l101 "There are no atoms"）
2. def2 基组名不能带连字符（def2-TZVP → def2TZVP，代码自动兼容）
3. formchk 归档源==目标 SameFileError（已处理跳过）

## 阶段 4：测试与验收 — ✅ 通过（2026-08-11）

| 用例 | 引擎 | 结果 |
|---|---|---|
| 乙醇 RESP（B3LYP/def2TZVP, pop=MK） | gaussian | sum=-0.000000；O=-0.6644、OH-H=+0.3924、C(CH2)=+0.5064、C(CH3)=-0.2583；校验全过 |
| 乙醇 RESP2（gas + water PCM, δ=0.5） | gaussian | sum=0；双单点 + RESP2 拟合；校验全过 |
| 苯 RESP（12 原子） | gaussian | 6C≈-0.109±0.002、6H≈+0.109±0.001（近似等价，见下）；校验全过 |
| QUICK 旧路径回归（乙醇） | quick | sum=0；9777 ESP 点两阶段；校验全过（零破坏） |
| 冒烟测试（甲烷/乙醇/混合） | — | 3/3 通过 |

**注意**：Multiwfn 标准 RESP 流程**不自动施加原子等价约束**（日志 "No atom
equivalence constraint is imposed"），苯电荷为近似等价（差异 <0.005 e，物理合理）。
严格等价需 7→18→5 eqvcons（可选增强，见 dev_notes 7）。

## 阶段 5：文档与归档 — ✅ 通过（2026-08-11）

- README：支持矩阵更新（RESP/RESP2 完成）、场景 2 重写（Gaussian 主路径 + qm 段）、FAQ 更新
- dev_notes：需求矩阵更新、4.1 Gaussian 主路径章节（含三坑 + 验证数据）、后续方向更新
- 提交：`b7cf6d4`（+728 行，11 文件）；代码清理 `b7cf6d4` 后 pyflakes 零告警（2026-08-11 整理）
