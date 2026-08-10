# smi2data 环境搭建文档（阶段 0）

**日期**：2026-08-09
**执行**：OpsAgent
**状态**：✅ 完成，全部验收通过

---

## 1. 环境信息

| 项目 | 值 |
|------|-----|
| 系统 | WSL2 Ubuntu 26.04 LTS（内核 6.18.33.2-microsoft-standard-WSL2） |
| 用户 | yryd（HOME=/home/yryd） |
| 硬件 | 8 核 / 11 GiB 内存 / 1 TB 磁盘 |
| Miniconda | `/home/yryd/packages/miniconda3`（conda 26.5.3） |
| 专用环境 | `/home/yryd/packages/miniconda3/envs/smi2data` |
| 环境大小 | 2.4 GB |
| Python | 3.12 |

> 所有软件均安装在 `/home/yryd/packages` 路径下（符合要求）。
> 原 `~/miniconda3`（含旧 md_env）保留未动，`.bashrc` 的 conda init 已改为指向新路径（备份：`~/.bashrc.bak.20260809`）。

---

## 2. 安装步骤（可复现）

### 2.1 安装 Miniconda 到 packages

```bash
# 安装脚本已在 /home/yryd/packages/Miniconda3-latest-Linux-x86_64.sh（197 MB）
chmod +x /home/yryd/packages/Miniconda3-latest-Linux-x86_64.sh
bash /home/yryd/packages/Miniconda3-latest-Linux-x86_64.sh -b -p /home/yryd/packages/miniconda3 -u
```

### 2.2 创建专用环境

```bash
/home/yryd/packages/miniconda3/bin/conda create -y \
  -p /home/yryd/packages/miniconda3/envs/smi2data \
  -c conda-forge \
  python=3.12 ambertools=26 rdkit parmed numpy pandas
```

### 2.3 更新 shell 配置

```bash
/home/yryd/packages/miniconda3/bin/conda init bash
```

### 2.4 激活环境

```bash
conda activate /home/yryd/packages/miniconda3/envs/smi2data
# 或（若 conda 未初始化）
source /home/yryd/packages/miniconda3/etc/profile.d/conda.sh
conda activate /home/yryd/packages/miniconda3/envs/smi2data
```

---

## 3. 安装清单（版本）

| 软件 | 版本 | 用途 | 来源 |
|------|------|------|------|
| antechamber | 26.0 | 原子类型(gaff2) + 电荷(AM1-BCC/ABCG2/RESP) | ambertools |
| parmchk2 | 26.0 | 缺失参数补齐 → frcmod | ambertools |
| tleap | 26.0 | 分子拓扑×N 拷贝 + 坐标绑定 → prmtop | ambertools |
| sqm | 26.0 | 半经验 QM（AM1-BCC 电荷） | ambertools |
| quick | 26.0 | 开源 QM（RESP 电荷） | ambertools |
| packmol | **21.0.1** | 装盒 | ambertools 自带 |
| parmed | 26.0 | prmtop → LAMMPS data | ambertools |
| rdkit | 2026.03.1 | SMILES → 3D | conda-forge |
| numpy | 2.4.6 | 数值计算 | conda-forge |
| pandas | 3.0.5 | 数据处理 | conda-forge |
| python | 3.12 | 运行环境 | conda-forge |

---

## 4. 关键经验与坑

### 4.1 ⚠️ 不要单独安装 packmol

AmberTools 26 的 conda-forge 元数据**自带 packmol**，并用 `packmol ==9999999999` 约束禁止安装外部 packmol 包。
直接 `conda install ambertools packmol` 会报：
```
LibMambaUnsatisfiableError: package ambertools-26.0 has constraint packmol 9999999999
```
**解决**：只装 `ambertools=26`，packmol 自动随附（实测为 21.0.1）。

### 4.2 ⚠️ WSL 后台任务保活

`wsl -e bash -c "nohup ... &"` 会在 wsl.exe 退出后把子进程杀掉（WSL2 空闲 VM 自动关闭）。
**解决**：Windows 侧用 `start /b wsl -e bash -c "... & wait"` 保持会话存活直到任务完成。

### 4.3 ⚠️ 网络断线恢复

安装中途网络中断会报 `CondaHTTPError: HTTP 000 CONNECTION FAILED` 并退出。
**解决**：网络恢复后重跑同一条 `conda create` 命令，已下载的包保留在缓存中，从断点继续。

---

## 5. 验收结果（阶段 0 全部通过）

### 5.1 命令行工具
```
antechamber: OK   parmchk2: OK   tleap: OK   packmol: OK
sqm: OK           quick: OK      parmed: OK
```

### 5.2 Python 库
```
rdkit 2026.03.1 | numpy 2.4.6 | pandas 3.0.5 | parmed: OK
```

### 5.3 ✅ ABCG2 实测（关键验证）

计划中担心 AmberTools 26 可能未含 ABCG2，实测**完全支持**：

```bash
antechamber -i test_ethanol.sdf -fi sdf -o test_abcg2.mol2 -fo mol2 -c abcg2 -at gaff2 -nc 0
```

- ✅ 无 "unknown charge type" 报错
- ✅ `-at` 选项列表明确包含 `abcg2`（gaff, gaff2, amber, bcc, abcg2, sybyl）
- ✅ 数据文件已随包安装：
  - `$CONDA_PREFIX/dat/antechamber/BCCPARM_ABCG2.DAT`
  - `$CONDA_PREFIX/dat/antechamber/ATOMTYPE_ABCG2.DEF`
- ✅ 乙醇 mol2 生成成功，GAFF2 类型正确（c3/oh/hc/h1/ho）
- ✅ **电荷守恒**：ABCG2 电荷总和 = 0.000001 ≈ 0（乙醇净电荷 0，容差 0.01 通过）

> **结论**：阶段 5（ABCG2 分支）无需补拷数据文件，可直接使用 `-c abcg2`。

### 5.4 测试文件
```
/home/yryd/packages/test_ethanol.sdf      # 乙醇 3D SDF（RDKit 生成）
/home/yryd/packages/test_abcg2.mol2       # ABCG2 电荷 mol2
```

---

## 6. 下一步（阶段 1：MVP 单分子 GAFF2 + AM1-BCC → data.lmp）

环境已就绪，可直接开始：
1. 编写 `smi2data` CLI 骨架（yaml 解析）
2. RDKit SMILES→SDF → antechamber(bcc,gaff2) → mol2 → parmchk2 → frcmod
3. ParmEd → data.lmp
4. 验收：苯/乙醇各出 data.lmp；LAMMPS 能读入
