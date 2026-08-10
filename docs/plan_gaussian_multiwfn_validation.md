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

## 阶段 2：安装 Multiwfn — 待执行

## 阶段 3：代码集成 — 待执行

## 阶段 4：测试与验收 — 待执行

## 阶段 5：文档与归档 — 待执行
