# Smooth4PC T73 — 条件 skein-lasagna 阻碍

[English](README.md)

本仓库支持对 trace-73 Cappell--Shaneson **同伦** \(4\)-球面
\[
A=\begin{pmatrix}0&269&1240\\0&41&189\\1&0&32\end{pmatrix}
\]
（Iwaki 标准形 \(X_{41,189,73}\)）的一份**条件性** skein-lasagna 阻碍论证。

**不断言光滑四维 Poincaré 猜想的反例。** 控制性文稿为
[`paper/spc4-t73-candidate/main.tex`](paper/spc4-t73-candidate/main.tex)
（中文对照：`main-zh.tex`；标题：《trace-73 Cappell--Shaneson 球面的条件
skein-lasagna 阻碍》）。

## 已证与开放

对显式 **Johnson 生成元** 柄表示，正文陈述 skein-lasagna 在 quantum degree
\(494\) 处比较所需的几何输入 **P0、C、S、P3**，以及识别
\(X_J\cong\Sigma_A^0\)。**P0、C、S 目前开放**：已提交的证书只检查这些输入的
有限组合模型（控制股、局部模型 movie、模型球面），不是实际的
Cappell--Shaneson 几何。论文第 3 节的状态表是本仓库的断言边界。

精确有限计算给出非零整数 \(2624\)，它是冻结 Burau 计算在 collar 端点约定下的
divided cubic；把它与实际几何类的 MWW divided cubic 等同是开放输入 C 的一部分。
Artin--Magnus 证书与 pure-braid Andreadakis 定理确立公开 braid 词的三阶性质。

Lean 开发将**抽象商论证**形式化：若给定将 MWW 商与四柄运输组装起来的接口数据
（`ExternalGeometry`），则非零次数-\(494\) 类**将**阻碍与 \(S^4\) 的微分同胚。
这些几何接口**未**在 Lean 中构造。

| 层次 | 状态 |
| --- | --- |
| 有限代数（\(2624\)、次数 \(494\)、\(\det A=\det(A-I)=1\)） | Lean 已检 |
| 抽象条件蕴含 | Lean 已检 |
| Johnson P0 / C / S（实际几何） | **开放**（仅有限模型证书） |
| P3 计算部分（E12：\(S^4\) 上空 link 次数 494） | 已检；E11/E13 依赖 P0 |
| Lean 中 `ExternalGeometry` 实例 | **开放** |

Lean 边界见
[`Smooth4PC/T73External.lean`](Smooth4PC/T73External.lean)。前提审计：

```text
python3 scripts/audit_t73_premises.py --check
```

预期摘要：`P0/C/S/P3=PASS`，`OVERALL=OPEN`，`COUNTEREXAMPLE=False`。
其中 `PASS` 是证书内部值（有限模型可重放）；论文状态表将 P0、C、S 记为开放。

**勘误（2026 年 9 月 2 日）。** 较早草稿混用两套 endpoint 索引表，报告了
\(-59072\)。统一到 braid 词所用 collar 表后，精确值为 \(+2624\)（仍非零）。

## 审阅 PDF

- 英文：[`output/pdf/spc4-t73-candidate.pdf`](output/pdf/spc4-t73-candidate.pdf)
- 中文：[`output/pdf/spc4-t73-candidate-zh.pdf`](output/pdf/spc4-t73-candidate-zh.pdf)

文稿与编译说明：
[`paper/spc4-t73-candidate/README.md`](paper/spc4-t73-candidate/README.md)。

## 从哪里开始

1. 阅读论文摘要与第 3 节（精确陈述）：
   [`paper/spc4-t73-candidate/main.tex`](paper/spc4-t73-candidate/main.tex)
   或上方英文/中文 PDF。
2. 按 [`REPRODUCING.md`](REPRODUCING.md) 从全新 clone 编译、检查公理报告并重算检测量。
3. 按下方命令重放有限检测量与 Johnson P0/C/S/P3 证书（完整清单见 `REPRODUCING.md`）。
4. 复核边界图见
   [`docs/INDEPENDENT_REVIEW.md`](docs/INDEPENDENT_REVIEW.md)。

## 重算 / 重放计算脚本

在仓库根目录使用 **Python 3.10+**。Windows 上可用 `python` 代替 `python3`。
每个 `--check` 会在内存中再生证书，并与 `audit/` 下已提交的 JSON 比对。

### 有限检测量（\(D_3=2624\)）

```text
python -I -B scripts/recompute_t73_delta3.py --check
```

预期：`DELTA3_ETA_T1=2624`、`DELTA3_XI=0`、`VERIFY=PASS`。

### Johnson P0 → C → S → P3（推荐顺序）

```text
# P0（约 1–2 分钟）：AR 桥、消去、几何辫
python -B scripts/certify_t73_p0_johnson.py --check

# 可选：写出并核验严格 P0 重构输入（默认未提交该 JSON）
python -B scripts/build_t73_p0_reconstruction_input.py --write
python -B scripts/reconstruct_t73_p0.py audit/t73_p0_reconstruction_input.json

# C：product 矩形、比较支撑、汇总 witness
python -B scripts/certify_t73_c1_cut_link.py --check
python -B scripts/certify_t73_c2_comparison.py --check
python -B scripts/generate_t73_c_comparison_witness.py --check

# S：反向 belt 球面与相对移动账本
python -B scripts/certify_t73_s_standard_spheres.py --check
python -B scripts/certify_t73_s_relative_moves.py --check

# P3：四柄图景、标准 S^4 次数 494、CS 识别
python -B scripts/certify_t73_p3_four_handle.py --check
python -B scripts/certify_t73_e12_s4.py --check
python -B scripts/certify_t73_e13_close.py --check
python -B scripts/certify_t73_e13_identification.py --check

# 前提汇总（须保持 OVERALL=OPEN / 非反例）
python -B scripts/audit_t73_premises.py --check
python -B scripts/check_t73_claim_boundary.py
```

| 脚本 | 作用 | 预期 |
| --- | --- | --- |
| `certify_t73_p0_johnson.py` | P0 Johnson 替换 | `T73_P0_JOHNSON_CERTIFICATE=PASS` |
| `reconstruct_t73_p0.py` | 严格 PL collar 对公开词 | `P0_RECONSTRUCTION=PASS`，`B44_LENGTH=11340` |
| `certify_t73_c1_cut_link.py` | 44 条带 + 227 剩余 \(z\) | `RECTANGLES=44`，`LEFTOVER_Z_CIRCLES=227` |
| `certify_t73_c2_comparison.py` | C2 支撑不相交 / \(H\) movies | `T73_C2_COMPARISON=PASS` |
| `generate_t73_c_comparison_witness.py` | C 账本绑定 P0/C1/C2 | `C_STATUS=PASS` |
| `certify_t73_s_*.py` | S 球面系与移动 | `PASS` |
| `certify_t73_p3_four_handle.py` | \(X_J\) 四柄层 | `E11`/`E12` PASS；`E13=PARTIAL` 属设计 |
| `certify_t73_e12_s4.py` | 标准 \(S^4\) 上空链次数 \(494\) | `S4_DEGREE_494_ZERO=True` |
| `certify_t73_e13_*.py` | \(X_J\cong\Sigma_A^0\) 管道 | `IDENTIFIED_WITH_SIGMA=True` |
| `audit_t73_premises.py` | 汇总状态 | `P0/C/S/P3=PASS`，`COUNTEREXAMPLE=False` |
| `check_t73_claim_boundary.py` | 论文保持条件性主张 | `T73_CLAIM_BOUNDARY=OPEN_GEOMETRY` |

**说明。**

- 证书经 SHA 链式绑定：C 绑定 P0，S 绑定 P0+C，P3 绑定 P0+C+S。改上游而不再生下游会导致 `--check` 失败。
- `certify_t73_p3_four_handle.py` 报 `E13=PARTIAL` 是预期行为；完整 \(\Sigma_A^0\) 识别在 `e13_*` 脚本中。
- C1/C2 核验的是正文所用的**组合 / PL 模型**义务，本身并不重证 MWW/BPW 范畴比较。
- Lean 编译（`tests/test_t73_minimal_formalization.py`）另计，较慢（约 5–10 分钟）；见 [`REPRODUCING.md`](REPRODUCING.md)。

## 复现约定

- Lean 工具链：`leanprover/lean4:v4.32.1`
- mathlib 修订：`520045ab14e26149ee970e2e617ca04b09bde5d6`
- Python：3.10 或更高
- 预期公理报告数：`38`
- 允许公理：`propext`、`Classical.choice`、`Quot.sound`
- `sorryAx`：无
- 预期检测量：`2624`

已提交的 `lake-manifest.json` 锁定 Lean 依赖。构建产物与本地依赖副本不在源码约定内。

## 范围

这是面向**条件性**阻碍的公开核验包，不是同行接受声明，也不是已证反例。
最有用的否定性审阅应针对 Johnson 几何绑定与尚未组装的 Lean
`ExternalGeometry`，而非已检整算术。

## 为何先在 GitHub 公开

GitHub 作为首个公开渠道，便于获取、速度与复现，不能替代学术审阅。既有
arXiv 投稿史在计算机科学方向；数学类目背书可能暂不具备。完整论证、Lean
源码、证书与重放说明均可在此检查。经实质审阅后，计划提交常规预印本。

## 许可证

本仓库以 [MIT License](LICENSE) 发布。
