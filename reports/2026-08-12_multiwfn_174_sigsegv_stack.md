# 检查报告：Multiwfn RESP2 大分子失败（退出码 174）诊断与修复

- **日期**：2026-08-12
- **任务**：MPC / SMBA（60/58 原子，RESP2，B3LYP/def2TZVP + PCM water）远程运行失败
- **现象**：`[warn] Multiwfn 退出码 174（可能仍已产出文件）` → `RuntimeError: Multiwfn RESP 未生成电荷文件: {name}_gas.chg`，MPC 与 SMBA 两个任务同点位失败；小分子对照（mpc-cut/smba-cut，32 原子）正常。

## 1. 结论（TL;DR）

- **不是**并发问题、**不是** g16 问题、**不是**下载/文件损坏问题。
- 根因：**Multiwfn 读大分子 fch（约 ≥1000 基函数，含 F 轨道）时，在"生成密度矩阵"阶段触发段错误（`forrtl: severe (174): SIGSEGV`），退出码 174，未产出 .chg**。
- 直接原因：Linux 默认进程栈上限 8 MB（`ulimit -s 8192`），Multiwfn 该步骤栈需求超限。Multiwfn 手册 §2.1.2 与官方论坛（id=207）明确：`ulimit -s unlimited` 解决。
- 修复：getlmp `_run_multiwfn` 在子进程内自动解除栈限制（`preexec_fn` + `RLIMIT_STACK=INFINITY`）+ 设 `OMP/KMP_STACKSIZE=1G`。**用户无需手动改环境**。

## 2. 证据链

| 检查项 | 结果 |
|---|---|
| g16 日志（mpc_gas/mpc_solv） | 均 `Normal termination`（11:25:56 / 11:32:03），C16H34NO6PS2，HF=-2274.62 |
| .chk 完整性 | WSL stat：mpc_gas.chk=64,397,312 B（log 报 Chk=62 MB ✓），smba 同理 ✓ |
| .fch 完整性 | 43,024,892 B，531,295 行，`Gaussian Version` 结尾；与 formchk 重新生成**字节级一致**（cmp IDENTICAL）✓ |
| 本地复现 | 同一 Multiwfn 2026.7.15 二进制 + 同一 fch → **稳定复现退出码 174**（与远程无关）|
| 崩溃点 | `Generating density matrix based on SCF orbitals...` 后 SIGSEGV |
| 小分子对照 | mpc-cut(548 基函数)/smba-cut：EXIT=0，正常产出 .chg ✓ |
| 排除项 | OMP_STACKSIZE=512M 无效（只管 OpenMP 线程栈）；OMP_NUM_THREADS=1 仍崩（非并行区问题）|

## 3. 修复验证（补丁后，不手动 ulimit）

shell `ulimit -s` 保持 8192，仅靠 `resp_from_fch` 内部 preexec_fn 解除限制：

| 文件 | 原子数 | 电荷数 | Σq | 结果 |
|---|---|---|---|---|
| MPC gas | 60 | 60 | 0.000000 | ✅ .chg 产出 |
| MPC solv | 60 | 60 | -0.000000 | ✅ |
| SMBA gas | 58 | 58 | -0.000000 | ✅ |
| SMBA solv | 58 | 58 | -0.000000 | ✅ |

RESP 两阶段拟合均收敛（Stage1 6 iter / Stage2 7 iter，含自动等价约束），Multiwfn 日志 `Result have been saved to {base}.chg`。

## 4. 代码改动

- `src/multiwfn.py`：
  - `_run_multiwfn` 增加 `preexec_fn=_unlimited_stack`（`resource.setrlimit(RLIMIT_STACK, INFINITY)`）
  - 子进程 env 设 `OMP_STACKSIZE=1G`、`KMP_STACKSIZE=1G`（手册 §2.1.2）
  - 模块 docstring + 函数注释记录该坑
- `docs/dev_notes.md`：§4.1 "Gaussian 链路实测三个坑" → 补第 4 坑
- commit：`148f48e`（main）

## 5. 远程机器处理（用户侧）

1. **推荐**：同步代码补丁（本仓库已提交，`git pull` 或按第 4 节手工改 `src/multiwfn.py`）。
2. **临时应急（不改代码）**：`start.sh` 里 `python main.py` 前加一行 `ulimit -s unlimited`。
3. 重跑 `MPC` / `SMBA` 两个任务即可（g16 单点结果可复用：chk/fch 完整，Multiwfn 会重新拟合）。
4. 备注：大分子 RESP 拟合耗时增加（ESP 网格计算 ~112 s/分子 + 拟合），属正常。

## 6. 遗留/建议

- Multiwfn 2026.7.15 该崩溃为环境（栈限制）触发，非版本 bug；不同机器默认 `ulimit -s` 可能不同，代码层修复已覆盖。
- 后续可在 `install.md` 环境配置里补充 `ulimit -s unlimited`（与论坛建议一致），双保险。
