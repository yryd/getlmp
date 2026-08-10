#!/bin/bash
# Gaussian 16 环境变量（供代码 subprocess 复用；source 后 g16 可用）
# 用法：source scripts/g16env.sh   （在项目根目录下执行）
#
# 注意：安装目录为 /home/yryd/packages/soft/g16/g16（tar 包顶层含 g16/ 目录，
# 解压到 soft/g16/ 产生嵌套；未做目录上移，g16root 直接指向实际安装目录）。
export g16root=/home/yryd/packages/soft/g16/g16
export GAUSS_EXEDIR=$g16root/bsd:$g16root/utility:$g16root
export GAUSS_SCRDIR=/home/yryd/packages/soft/scratch
export LD_LIBRARY_PATH=$g16root/bsd:$g16root
export PATH=$g16root:$PATH
