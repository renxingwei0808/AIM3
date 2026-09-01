# 任务：AIM3 仓库服务器批次 1（Tier B/C）—— 由本地 Claude Code 通过 SSH 驱动远端服务器

你在本地 WSL 运行，仓库 `renxingwei0808/AIM3` 的 `main` 已包含 Tier-A 批次 1 的全部结果（`results/REPORT_tierA_batch1.md`）。请先读 `CLAUDE.md`、`docs/ROADMAP.md`（含文末 Status update）、`server/README.md`、`server/SIZING.md`、`results/REPORT_tierA_batch1.md` 第 4、6、7 节，以及本批次涉及的实验 README。本批次的目的：把 Tier A 上"[VERIFIED n=11,14；CONJ 一般情形]"的 Gate 1 结论推进到跨多个 n 的验证或证伪，并首次测量真实的 F4 求解次数。

## 0. 工作方式（硬性要求）

1. 分支 **`claude/server-batch-1`**，从 `main` 新建。本批次只写以下位置，避免与并行的网页端分支冲突：`results/server/**`、`results/E*_server*.csv`（新文件）、`results/REPORT_server_batch1.md`、各实验 README 中新增的 "## Server log" 小节、`experiments/*/magma|sage/*` 与你新增的后端脚本、`server/`。**不改** `docs/`、`paper/`、`legacy/`、已有的 Python 实验脚本主体（允许增加 `--backend` 选项，见任务 2）。
2. 服务器上所有长任务在 `tmux` 中用 `nohup … > log 2>&1 &` 启动，轮询日志，不阻塞会话；一次只跑一个内存 > 50 GB 的作业，`ulimit -v` 按 `server/README.md` 设置；每个作业记录 `host, commit, seconds, peak_RSS`（用 `/usr/bin/time -v`），超时/超内存只报下界并标 `censored=1`。
3. 任何 Magma/Sage 脚本在用于 n > 13 之前，必须先在 n ≤ 13 上与 `results/E03_sbox_degree_fall.csv`、`results/E12_first_layer_ideal.csv`、`results/E06_kernel_model_falls.csv` 中已有的行逐项对齐（秩、syzygy、降次全部一致），把对齐结果写进 Server log。
4. 只把摘要 CSV/日志同步回仓库；矩阵、中间基不入库。commit 按任务分开；最后 push 分支并给出 PR 描述。
5. 状态标签规则同 `CLAUDE.md`；与 Tier-A 结论矛盾的结果照实报告并单独成节。

## 1. 环境（第一天）

在服务器上：`git clone` 仓库到工作目录；确认 `magma` 可用并记录版本与 `SetNthreads` 上限；安装/确认 Sage ≥ 10（优先 conda-forge 或已有安装）、msolve（源码编译，`-t` 线程）、GNU parallel；跑 `bash scripts/smoke_test.sh` 与 `python3 tests/run_tests.py`（必须 6/6）。把软件版本写入 `server/ENVIRONMENT.md`（新文件）。

## 2. 快速秩后端（任务 3–5 的前提）

Tier-A 的 `src/aim3/boolmac.py` 是纯 Python 流式消元，在 n=13、D=4（约 9×10⁴ 行 × 9×10⁴ 列）上需要数天。请增加一个不改变语义的后端：
- 新建 `src/aim3/boolmac_export.py`：把 Macaulay 矩阵按 `boolmac.py` 的列序（D 次单项式在最高位）导出为稀疏行文件；
- 新建 `experiments/E04_full_system_hilbert/sage/rank_backend.sage`：读入后，用 Sage 的 GF(2) 稠密矩阵（底层 M4RI）计算 `rank(M)` 与 `rank(M restricted to degree-D columns)`，输出 rank、syz、`dim(span∩deg≤D−1) = rank(M) − rank(top block)`；
- 在 `full_system_falls.py`、`first_layer_ideal.py`、`kernel_bilinear_toy.py`、`sbox_degree_fall.py` 增加 `--backend sage`（默认仍为 python），调用上述导出 + Sage；
- 验证：对 (E03: n=11 e=3; n=13 e=7)、(E12: n=11 --D4)、(E06: n=7 D=4) 两种后端输出完全一致后，才允许在更大 n 上使用。

## 3. Gate 1 最终判定（最高价值）

`full_system_falls.py --backend sage --per-block`：
- ℓ=2（三个 S 盒，AIM3-I/III 形状）：n=13、15、17，seeds 1,2,3，系统 S_quad / P1only / random 控制；
- ℓ=1：n=13、15；ℓ=3（四个指数，AIM3-V 形状）：n=11、13；
- 每次运行报告饱和度 `dim_low/C(nv,≤3)`；只有饱和度 < 0.9 的运行进入判定。判定量是 `cross_excess = falls(all) − Σ falls(单 S 盒嵌入同变量环)`；同时记录 `falls(all)` 与 `(ℓ+1)·2n` 供对照。
- 结论按 `COMPLEXITY_ACCOUNTING §4`：若在 ≥ 4 个 n 上 cross_excess = 0（允许退化指数造成的 ≤ 1 残差），Gate 1 记为 "直和 [VERIFIED n=…]"；若 cross_excess 随 n 增长，记为 "存在跨 S 盒降次 [VERIFIED]"，并给出增长律与涉及的 S 盒对（用 `--per-block` 的两两组合定位）。这一节在报告中单独成节，写明对 T1/T2 优先级的含义。

## 4. 降次规则在拟合范围外的检验（证伪优先）

`bash experiments/E03_syzygy_and_degree_fall/sweep.sh 17 18 19 20`（`--backend sage`），然后 `summarise.py --md`。规则 `predicted_falls_D4` 是对 n ≤ 16 的拟合：报告 n=17..20 全部合法 e 的 mismatch 数；任何 mismatch 都要给出 (n, e, 实测, 预测, 标签)。若无 mismatch，把 `paper/claims_ledger.md` 中该项的建议标签写在报告里（仍为 [CONJ]，但范围扩到 n ≤ 20）——不要自己改 ledger。另外跑 n=14、16 的 D=5（`SIZING.md` 估计 8–20 GB），报告 D=5 是否出现新的降次到 ≤4 次，作为 T1.2 的补充数据。

## 5. E12 与 E06 的延伸

- `first_layer_ideal.py 13 --D4 --backend sage` 与 `14 --D4`；`cross_relations.py 16 4` 后运行 `explain_n16.py 16 --deg 4`（若脚本不支持 deg 参数，按其 README 的方法先扣除单变量像判据 Tr(t_i^{−r_i}) 再报告跨 S 盒四次关系维数）。
- `kernel_bilinear_toy.py n 4 --backend sage`，n = 11, 13, 15（含随机双线性对照）：报告 D=3 新二次式是否恒为 #eqs = 6n，D=4 新三次式的数量与 n 的关系（拟合前先看饱和度）。

## 6. 真实求解次数：F4 步次数日志（回答 E04 Q2）

这是饱和不会破坏的观测量，也是本批次真正新的东西。
- 用 `legacy/aim3_macaulay.py` 的 `build_aim3` 为 n = 7, 9, 11, 13, 15、ℓ = 1, 2、seeds 1,2 生成 S_quad / P1only / random 控制，导出为 msolve 输入格式与 Magma `BooleanPolynomialRing` 格式（新增导出脚本 `experiments/E04_full_system_hilbert/export_systems.py`，输出到 `results/server/<date>_E04/systems/`）。
- msolve：`msolve -v 2 -t 32 -f sys.ms -o out.ms`，解析日志得到每一步的次数与矩阵尺寸，取最大步次数为 `d_solv`；Magma：`SetVerbose("Faugere",2); GroebnerBasis(...)`（`experiments/E04_full_system_hilbert/magma/aim3_f4_degrees.m` 是骨架，补全并在 n=7 与 msolve 交叉验证）。
- 同样对 E06 核模型（`kernel_bilinear_toy.py` 增加 `--export`）和 E10 的块序（`elim` 序：pt 块在前 / t 块在前；核模型 U>W / W>U，Magma）跑一遍。
- 输出 `results/E04_f4_degrees.csv`：n, ℓ, exps, seed, system, order, tool, d_solv, max_rows, time, peak_RSS。报告中给出 d_solv vs n 的表（不做外推；只有 ≥ 4 个 n 才允许写增长趋势），并与 B2 假定的 d_reg 对照。

## 7. E05：降变量系统的 resultant（Magma，Tier C）

按 `experiments/E05_two_valued_inverse/README.md` 的 "Input format for the server"：用 `reduced_system_toy.py --export`（若无此选项则增加，输出两条稀疏方程为 Magma 可读格式）在 n = 8..20（seed 1, 正确分支）生成实例，Magma 中 `Resultant(Eq1, Eq2, t1)`，记录：耗时、内存、结果多项式次数与项数、根数与真密钥是否为根。拟合 `log2(time)` 与 `log2(deg)` 对 n 的斜率；ROADMAP §3.3 的终止判据是斜率 ≥ 2。若 Magma 在某个 n 超时/超内存，报下界并停止上探。若时间允许，测试 Magma `Factorization` 或 `Decomposition` 是否让结果多项式分解（T3.3 的线索）。

## 8. 交付

`results/REPORT_server_batch1.md`：① 环境与后端对齐记录；② Gate 1 判定（独立成节）；③ 降次规则 n=17..20 的检验；④ E12/E06 延伸；⑤ d_solv 表；⑥ E05 resultant 曲线与斜率；⑦ 与 Tier-A 结论矛盾或无法解释的现象；⑧ 未完成项；⑨ 对 `paper/claims_ledger.md` 的建议行（不直接改）。push 分支 `claude/server-batch-1`，给出 PR 描述（改动文件、数字总表、待人工决策事项）。

