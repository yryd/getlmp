# getlmp 环境安装与配置指南

> 本文档是 getlmp 运行环境的**唯一权威安装说明**。代码本身不做任何路径探测
> （不写死候选目录、不自动 source 脚本），只按下列**约定**工作。安装一次后，
> 在任何目录运行 `python main.py input.yaml` 都能找到所需命令。

## 0. 环境约定（代码如何找工具）

| 依赖 | 约定 | 代码探测方式 |
|---|---|---|
| Python 库（rdkit/parmed/yaml/numpy） | 安装在 conda 环境内 | 直接 import |
| antechamber / parmchk2 / tleap / sqm / resp / respgen / packmol / quick | conda 环境（ambertools） | `shutil.which()`（PATH） |
| g16 / formchk | 任意目录，**必须加入 PATH**（.bashrc） | `shutil.which()`（PATH） |
| Multiwfn_noGUI | 任意目录，**必须加入 PATH**（.bashrc） | `shutil.which()`（PATH） |
| `QUICK_BASIS` 环境变量 | 指向含 `6-31G.BAS` 的目录（.bashrc） | 读环境变量并校验 |

**要点**：
- 代码里**没有**任何硬编码路径；所有外部命令一律从 PATH 找。
- 环境变量（`GAUSS_*`、`QUICK_BASIS`）由 `.bashrc` 提供，子进程直接继承。
- yaml 配置仍支持显式覆盖：`qm.g16root`、`qm.multiwfn_path`
  （一般无需配置，仅用于非常规安装位置）。
- 安装好之后用 `python check_env.py` 自检，缺什么一目了然。

## 1. 安装 conda 环境（RDKit / AmberTools / ParmEd）

```bash
conda create -n smi2data -c conda-forge python=3.12 \
    ambertools=26 rdkit parmed numpy pyyaml -y
conda activate smi2data
```

- `ambertools=26` 自带：`antechamber`、`parmchk2`、`tleap`、`sqm`、
  `resp`、`respgen`、`packmol`、`quick`。
- **坑**：不要单独安装 packmol —— AmberTools 26 自带 packmol 并用约束
  `packmol ==9999999999` 禁止外部安装，`conda install ambertools packmol`
  会报 UnsatisfiableError。
- 网络断线报 `CondaHTTPError` 后重跑同一条 `conda create` 命令即可续装。

## 2. 安装 Gaussian 16（RESP/RESP2 默认引擎 `qm.engine: gaussian`）

> 不需要 license 也能跑本地单点（G16 本地计算无需联网校验，详见 g16 包内安装说明）。
> 也可以不装 G16，直接用第 4 节的 QUICK 回退路径（但无 RESP2、无隐式溶剂）。

```bash
# 1) 解压到任意目录（示例 ~/packages/soft/g16/，顶层含 g16/ 子目录）
tar xjf /path/to/G16-C01-AVX.tbJ -C ~/packages/soft/g16/

# 2) 建 scratch 目录（Gaussian 运行时写临时文件）
mkdir -p ~/packages/soft/scratch

# 3) 把下面的环境变量写入 ~/.bashrc（注意改成你的实际路径）
```

**`~/.bashrc` 追加**（路径按你的安装位置修改）：

```bash
# ---------- Gaussian 16 ----------
export g16root=$HOME/packages/soft/g16/g16
export GAUSS_EXEDIR=$g16root/bsd:$g16root/utility:$g16root
export GAUSS_SCRDIR=$HOME/packages/soft/scratch
export LD_LIBRARY_PATH=$g16root/bsd:$g16root
export PATH=$g16root:$PATH
```

## 3. 安装 Multiwfn（RESP/RESP2 拟合 + ESP 可视化）

```bash
# 解压到任意目录（示例 ~/packages/soft/multiwfn/）
unzip /path/to/Multiwfn_2026.7.15_bin_Linux_noGUI.zip -d ~/packages/soft/multiwfn/

# ~/.bashrc 追加（noGUI 版可执行文件是 Multiwfn_noGUI）
export PATH=$HOME/packages/soft/multiwfn/Multiwfn_2026.7.15_bin_Linux_noGUI:$PATH
```

## 4. 配置 QUICK 回退路径（`qm.engine: quick`）

QUICK 随 ambertools 安装，但**基组数据目录**需要环境变量告知：

```bash
# ~/.bashrc 追加（$CONDA_PREFIX 指向 smi2data 环境）
export QUICK_BASIS=$CONDA_PREFIX/AmberTools/src/quick/basis
```

`QUICK_BASIS` 必须指向含 `6-31G.BAS` 的目录，代码会校验；缺失时给出友好报错。

## 5. 完整 `.bashrc` 模板（当前机器实测可用）

```bash
# ---------- conda ----------
source ~/packages/miniconda3/etc/profile.d/conda.sh
conda activate smi2data

# ---------- Gaussian 16 ----------
export g16root=$HOME/packages/soft/g16/g16
export GAUSS_EXEDIR=$g16root/bsd:$g16root/utility:$g16root
export GAUSS_SCRDIR=$HOME/packages/soft/scratch
export LD_LIBRARY_PATH=$g16root/bsd:$g16root
export PATH=$g16root:$PATH

# ---------- Multiwfn ----------
export PATH=$HOME/packages/soft/multiwfn/Multiwfn_2026.7.15_bin_Linux_noGUI:$PATH

# ---------- QUICK 基组 ----------
export QUICK_BASIS=$CONDA_PREFIX/AmberTools/src/quick/basis
```

改完执行 `source ~/.bashrc`（或重开终端）生效。

## 6. 自检

```bash
python check_env.py               # 全量
python check_env.py --quick       # 只查 QUICK 场景
python check_env.py --gaussian    # 只查 Gaussian 场景
```

输出 `✓`/`✗` 状态表；缺失项会给出修复提示。退出码 0=必需项就绪，1=有缺失。

## 7. 快速验证

```bash
# conda 命令（应在 PATH）
antechamber -h | head -1 && parmchk2 --help | head -1 && packmol < /dev/null | head -1

# Gaussian（应在 PATH，跑通一个水单点）
echo -e '%nprocshared=4\n%mem=2GB\n#p hf/6-31g*\n\ntest\n\n0 1\nO  0.0 0.0 0.0\nH  0.0 0.0 0.95\nH  0.95 0.0 0.0\n' > h2o.gjf
g16 h2o.gjf h2o.log && grep 'Normal termination' h2o.log

# Multiwfn（打印版本即成功）
Multiwfn_noGUI 2>&1 | head -5

# QUICK 基组
ls $QUICK_BASIS/6-31G.BAS
```

## 8. 常见问题

| 现象 | 原因与解决 |
|---|---|
| 报错 `PATH 中找不到 g16` | G16 未安装或 PATH 未配；按 §2 配置 .bashrc |
| 报错 `PATH 中找不到 Multiwfn_noGUI` | 按 §3 加入 PATH，或 yaml 配 `qm.multiwfn_path` |
| 报错 `QUICK_BASIS 未设置或无效` | 按 §4 配置；确认目录含 6-31G.BAS |
| g16 运行卡住/异常 | 确认 `GAUSS_SCRDIR` 目录存在且可写 |
| 非交互 shell 找不到 conda 命令 | 先 `source <conda>/etc/profile.d/conda.sh && conda activate smi2data` |

---

> 历史说明：早期版本曾用项目内 `scripts/g16env.sh` 封装 G16 环境变量并在代码里
> 注入子进程。2026-08-11 起废弃该方案——**安装方法统一收敛到本文档**，代码只认
> PATH 与 `.bashrc` 环境变量，删除 `scripts/` 目录与所有硬编码候选路径。
