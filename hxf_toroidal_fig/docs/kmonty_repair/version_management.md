# kmonty 双工作目录与实验管理

2026-10-06 按用户定义：历史 tools/RaptorP/kmonty 就是本项目的“原版/稳定版”。这是版本命名，不新增物理验收声明。

| 用途 | qposim分支/标签 | 工作目录 |
|---|---|---|
| 稳定实验源码 | codex/kmonty-stable / kmonty-stable-20261006 | /mnt/mydisk/qposim-kmonty-stable/hxf_toroidal_fig/postprocess/radiation/kmonty |
| κ开发源码 | codex/kmonty-pump-fix / kmonty-kappa-checkpoint-20261006 | /mnt/mydisk/qposim/hxf_toroidal_fig/postprocess/radiation/kmonty |
| 旧κ分支 | codex/kmonty-g0-g1 | 留在原66fbbd39作历史记录，不作为最新开发入口 |
| 实验文件 | 不在源码工作树内运行 | /mnt/mydisk/qposim-kmonty-experiments/stable/<实验名> 或 kappa/<实验名> |

两处源码由同一个qposim Git仓库管理，工作区、index、对象构建目录各自独立；Git提交对象和分支引用共享。不要在原工作目录来回checkout，也不要同时在两个工作目录checkout同一个分支。

稳定分支导入RaptorP提交538a0412f3324a2cd961cccf0e6ae5aabcd049c3的72个kmonty受跟踪文件；逐字节保留，包括现有src_var_kappa。旧tools/RaptorP/kmonty保持原样，仅作来源。开发分支将当前维护目录71个源码/配置/测试文件原样提交；既有证据仍保留原路径，绑定哈希见version_baseline.json。旧README中的历史科学状态不作为新的验收结论。

## 使用规则

1. 稳定版实验从稳定工作目录构建；κ修改只在pump-fix工作目录提交。基线标签不移动，之后的修改产生新提交/新版本标签。
2. 改编译宏、源码或架构后显式重新构建，使用新的独立OBJDIR。稳定版旧makefile的默认GSL相对路径不适用于新位置：两边均显式传入 `GSL_PREFIX=/mnt/mydisk/qposim/tools/raptorp-deps/prefix/usr`，并使用 `/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/comm_libs/mpi/bin` 工具链。
3. 运行前把确定的可执行文件复制进独立实验目录，不能链接到开发目录中会被覆盖的二进制。源文件、对象目录、输出目录和可写缓存不在两条路线之间共享。
4. 每个实验保存：源码commit、工作树是否dirty及必要diff、编译命令/有效宏、MPI/GSL版本、二进制SHA256、输入配置与快照路径、seed、MPI/OMP布局。运行中不要覆盖这份二进制或配置。
5. grid.out等只读模拟输入可以使用明确的只读来源路径/链接；Compton表、权重表、spec输出等可写文件留在各自实验目录。旧缓存只有身份完全匹配时才可显式复制，不能共用可写缓存。
6. 构建成功后先保留所需二进制，再按项目规则清理该次独立构建目录；先检查make clean作用范围，不能清理运行结果。
7. κ修复成熟后按具体提交审查、验证，再决定是否移入stable；不整条合并历史RaptorP G3工作树，不把κ物理改动自动加入stable。

查看版本（只读）：

```bash
git -C /mnt/mydisk/qposim worktree list
git -C /mnt/mydisk/qposim log -1 --oneline codex/kmonty-pump-fix
git -C /mnt/mydisk/qposim log -1 --oneline codex/kmonty-stable
git -C /mnt/mydisk/qposim diff codex/kmonty-stable..codex/kmonty-pump-fix -- hxf_toroidal_fig/postprocess/radiation/kmonty
```

本次只完成本地版本管理和哈希校验，没有重编译、重跑模拟或改变物理源码，也没有将稳定版命名解释为新生产验收。qposim未配置远程；本地worktree不是异机备份。将来配置私有远程后再明确推送分支和标签。
