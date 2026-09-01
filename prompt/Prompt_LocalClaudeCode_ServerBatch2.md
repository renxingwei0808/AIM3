# 任务：AIM3 仓库服务器批次 2（本地 Claude Code → 远端服务器）

批次 1 已合并到 `main`；`results/server/EXTRA_EXPERIMENTS.md` 与 `REPORT_server_batch1.md` §8 列出的仍在运行/未收集的作业先收集入库。请读 `REPORT_server_batch1.md` 全文、`REPORT_tierA_batch2.md` §2/§5、`docs/THEORY_NOTES.md`。分支 **`claude/server-batch-2`**。写入范围同批次 1（`results/server/**`、`results/E*_server*.csv`、`results/REPORT_server_batch2.md`、各实验 README 的 Server log、`experiments/*/magma|sage/*`、`server/`、E03/E04/E06/E12 脚本的向后兼容选项、`src/aim3/boolmac.py`）；**不动** `docs/`、`paper/`、`legacy/`、E07/E09（网页端分支正在写）。作业规则同批次 1（tmux/nohup、`/usr/bin/time -v`、内存上限、censored 标注、先小 n 对齐再放大）。

本批次回答三个问题：**(Q1) 校准的 solving-degree 律是什么；(Q2) S 盒理想在 4 次及以上是否还有二次理想之外的内容；(Q3) 核模型与 twisted-bilinear 结构能否被非 Boolean 的求解策略利用。**

## 0. 收尾（先做）
收集批次 1 仍在运行的作业（ℓ=3,n=15 三首层子集；E06 核模型块序 n=7；D=5 随机对照 n=16；ℓ=2 n=13/15 剩余 F4 行；Gate-1 的 seeds 2,3），入库并更新 `E04_server_gate1*.csv`、`E04_f4_degrees.csv` 等。跑批次 1 §8 提到的**解耦对照**：四个 S 盒放在不相交变量块上（无共享 pt），ℓ=3, n=11 子集格，确认 cross_excess=0，从而证明 n=11 的 1639 不是测量本身的问题。

## 1. 指数修正与校准律（Q1，最高优先）
1. 在 `full_system_falls.py`、`export_systems.py`、`kernel_bilinear_toy.py` 增加 `--generic-exps`：选取满足 `aim3.structures.degenerate_class(n,e)==[]` 的最小 ℓ+1 个合法 e（若不足则报错）。所有本批次的新实例都用它。
2. **d_solv 重跑**（Magma 为主，msolve `-g 2` 交叉）：S_quad / P1only / random，ℓ=2，n = 9, 11, 13, 15，seeds 1,2；ℓ=1，n = 11, 13, 15, 17（ℓ=1 变量少，n=17 可能可行；峰值内存 >150 GB 即放弃并记录）。
3. **混合猜测下的 δ**：对 ℓ=2 的 n=13, 15，固定 pt 的前 k 位（k=2,4,6，用真值），测 d_solv；同时用 `E02/table4_reproduce.dreg_boolean(3n−k, 9n)` 给出 B2 的预测。报告 δ(n,k) = d_reg^{B2} − d_solv 表。
4. **Hilbert 函数侧证**：对 ℓ=2, n=13 的 S_quad 在 D=5 计算 Macaulay 秩（`SIZING.md`：约 96 GB，单独作业），与布尔半正则级数在 D=5 的预测系数比较；D=5 的 cross_excess（直和是否在 D=5 仍成立，用 `gate1_pairs.py` 的子集格，只做单点与全集）。
5. **校准数字**：新增 `server/jobs/calibrated_estimate.py`：对 AIM3-I/III/V，用布尔级数求 d_reg(N−k, 9n)，按 δ∈{0,1} 两种情形取 d = d_reg − δ，成本 2^k·(Σ_{i≤d} C(N−k,i))^ω，ω∈{2, 2.81}，对 k 取最小；输出与 B1、B2、B3 的对照表。**只输出表，不在报告里做安全结论**——那是论文写作时的判断。

## 2. S 盒理想的高次结构（Q2）
1. `sbox_degree_fall.py` 增加 `--pointset-deg 4`：在 n=16、18 直接计算 dim I(V_e)_{≤4}（评估矩阵 65535×41449 / 262143×66711，M4RI），与 D=5、D=6 的 Macaulay 可达四次部分比较（D=6 在 n=16：行 48·C(32,≤4)≈2.0M，列 C(32,≤6)≈1.15M，约 290 GB，**超限则改用 n=14 的 D=6 或只报 D=5**）。报告"二次理想之外的四次生成元"数量；对 e=2 类同样做（已知有 2n 个三次式不可达）。
2. D=5 降次的分解：对 n=14,16（一般类与 {−3} 类），算 `span{z_k·c_j}`（2n 个新三次式 × 2n 个变量）在 degree-≤4 空间中的维数，与 falls_D5 相减得到"净新四次式"；验证批次 1 的 n(4n+7) 是否 = (2n)(2n) + 7n 的分解，以及 {−3} 类的 n(4n−6) 是否意味着三次式的倍式之间有 8n 个线性关系。
3. 若 (1) 发现理想外的四次生成元：输出它们的显式基（保存为 JSON），并在 ℓ=2, n=13 把它们加入 S_quad 后重测 d_solv（这是 Q2 通向攻击的唯一路径）。

## 3. 核模型的非 Boolean 求解策略（Q3）
1. **W-限制 XL**：在 `kernel_bilinear_toy.py` 增加 `--mult-block W|U`：Macaulay 乘子只取 W（或只取 U）块的单项式，测 D=2..5 的秩与降次，n=7,9,11；与全乘子对照。问题：只乘 W 是否在更低的 D 达到线性化（秩 = 单项式数 − 1）？
2. **MinRank/秩缺陷测试**：把第一层写成 N(W)·(U,1)ᵀ=0（N 为 2n×(n+1)，元素为 W 比特的仿射函数）。n=7,9：对真值 W 验证 rank N = n；对随机 W 统计 rank 分布；对 KS 建模（核向量的 Boolean 变量）用 msolve 求解，记录 d_solv 与时间，与 S_quad 对照。
3. **Sage/PolyBoRi 块序**（批次 1 的遗留）：用 Sage 的 `BooleanPolynomialRing(order='block')`（PolyBoRi 直接在块序中计算）对 n=7,9 的 S_quad 与核模型跑 GB，记录最大步次数（PolyBoRi 的 verbose 或按 `groebner_basis(deg_bound=…)` 逐次探测），比较 grevlex 与 U>W、W>U。若 PolyBoRi 不可用则说明并跳过。

## 4. T4 交叉验证（与网页端分支配合）
网页端分支将提交纯 Python 的 Ore 环原型与 `results/E07_theory_skew_degrees.csv`（n≤12）以及实例 JSON。本批次用 Sage `OrePolynomialRing(K, K.frobenius_endomorphism())` 在 n=8..24 上重算 lclm/gcrd 次数与系数关于 t₀ 的次数 D(n)（同样用求值-插值），先在网页端已有的 n 上对齐，再放大。**若网页端分支尚未合并，本任务顺延到其合并后，先做 1–3。**

## 5. E05 廉价消元（预算 4 小时）
按批次 1 §8 的两条建议：(a) 指数模 2ⁿ−1 后做 Kronecker 代换再 `Resultant`；(b) 利用 Eq1 的可分形式 W(t₀)+V(t₁) 用查表法消去 t₁。目标是取得 n=7..10 至少 3 个完成的尺寸，报告 log₂(time) 与 log₂(deg) 对 n 的斜率（≥3 点才报斜率）。

## 6. 交付
`results/REPORT_server_batch2.md`：收尾与解耦对照；Q1 的 d_solv 表（一般类指数）、δ(n,k) 表、D=5 Hilbert 侧证、校准数字表；Q2 的理想外生成元结论与 D=5 分解；Q3 的三项结果；T4 交叉验证（若做）；E05 斜率；矛盾与未完成项；ledger 建议行（不直接改）。push `claude/server-batch-2`，给出 PR 描述。