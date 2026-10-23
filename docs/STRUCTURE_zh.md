# 仓库结构与使用说明（中文）

本仓库是博士论文 *Institutional Constraints and Public Data Commercialization: Evidence from a Chinese Meteorological Media Enterprise*（Teesside University，白宇）的**复现材料包**，对应论文 4.8.4 节、4.9.4 节、附录 C 和附录 H.4 中的承诺：

> 完整分析代码（Jupyter/Python、Mplus 输入文件、requirements.txt）与一份**模拟数据集**公开发布于本 GitHub 仓库，并通过 Zenodo 存档获得持久 DOI；去标识化的个体层数据受附录 H.2 所述保密协议限制，不公开存放，按合理请求提供。

## 一、目录说明

| 目录 / 文件 | 内容 | 对应论文位置 |
|---|---|---|
| `README.md` | 英文总说明（审稿人、国际读者首先看到的文件） | — |
| `LICENSE` | 代码 MIT；模拟数据与文档 CC BY 4.0 | — |
| `CITATION.cff` | 引用元数据，GitHub 与 Zenodo 自动读取（作者、题名、年份、ORCID 可补） | — |
| `requirements.txt` | Python 3.11 依赖版本（与附录 C 开头一致） | 附录 C |
| `code/00_generate_synthetic_data.py` | 按论文公开参数（Table 5-1 均值/标准差/相关、Table 4-4 载荷、Table 5-3 路径、ICC）生成模拟面板 | Table 5-1、5-2、5-3、§4.9.3 |
| `code/01_data_preprocessing.py` | 数据质量检查：范围、缺失、重复、VIF、\|Z\|>3、Mahalanobis、Cook's D | 附录 C.1；§4.4.3、§5.4 |
| `code/02_measurement_checks.py` | α、CR、AVE、Fornell–Larcker、HTMT、第 1 波 EFA（特征值、KMO、Bartlett、平行分析）、ICC 与有效样本量、Harman 单因子、标记变量 | §4.5.7、§4.7.1、§4.9.3、§5.4、Table D-1 |
| `code/03_structural_paths_crosscheck.py` | 复合指标路径回归（员工聚类稳健标准误）、单元×年份双向固定效应 | 附录 C.2；Table 5-3、5-11、附录 D.2 |
| `code/04_bootstrap_mediation.py` | statsmodels Mediation（附录 C.3 原样）+ 按员工整体重抽样的聚类 bootstrap，同时报告 percentile 与 bias-corrected 区间 | Table 5-4、5-16；§4.6.1 |
| `code/05_iv_estimation.py` | 2SLS：IV1 邻省政策支持、IV2 观测站密度；第一阶段 F、Sargan J、Wooldridge 内生性检验 | 附录 C.4；Table 5-3 Panel B；§6.5.2 |
| `code/06_moderation_interaction.py` | H10 PR1×PR2 层次回归、简单斜率、Johnson–Neyman；H8/H9 按 DA/ET 分组的 a、b、c′ 与间接效应差异 bootstrap | Table 5-6、5-7、5-13；§5.5.2–5.5.3 |
| `code/07_robustness_checks.py` | 留一年、留一单元、替代测量、安慰剂置换、缩尾、剔除 Cook's D、分区域、按单元的 wild-cluster bootstrap | 附录 C.5；Table 5-5、5-8–5-12、D-3、D-4 |
| `code/08_mmdv_worksheet.py` | MMDV 工作表：CBV / MBV / IBV、等权与加权综合值、三情景敏感性（结果与附录 E 完全一致） | 附录 E；Table 5-14、5-15 |
| `code/09_prepare_mplus_data.py` | 导出 Mplus 用的无表头 `panel.dat` 与变量顺序文件 | — |
| `code/10_sem_latent_semopy.py` | **不用 Mplus 的潜变量 SEM**：用 semopy 拟合与 mplus/01–05 完全相同的五个模型，员工聚类 bootstrap 给出标准误、percentile/BC 区间、简单斜率、ET/DA 组间差异 | Table 5-2、5-3、5-4、5-6、5-7、5-13、5-16 |
| `code/11_prepare_real_panel.py` | 把作者本机的原始问卷/档案 Excel 转成编码手册格式（只在本机运行，输出被 .gitignore 排除） | — |
| `code/run_all.py` | 一键按顺序运行 00–09，日志写入 `outputs/run_log.txt` | — |
| `code/utils.py` | 公共函数（α、CR/AVE、HTMT、ICC、KMO、聚类 bootstrap、BC 区间） | — |
| `mplus/*.inp` | Mplus 8.8 输入：01 二阶 CFA；02 结构模型（MLR、TYPE=COMPLEX）；03 5000 次 bootstrap 中介（BCBOOTSTRAP 区间）；04 PR1×PR2 交互；05 ET/DA 多组模型 | §4.6.1；Table 5-2、5-3、5-4、5-6、5-7、5-13 |
| `mplus/README.md` | 各 .inp 用途、运行顺序、变量名缩写对照 | — |
| `data/synthetic/` | 模拟数据（1,800 人-波 × 90 变量；54 单元-年）、生成参数、达标日志 | — |
| `data/real/README.md` | 说明真实数据不在仓库内及申请方式（**不要**把真实数据放进来） | 附录 H.2、H.4 |
| `docs/CODEBOOK.md` | 逐变量编码手册：题项代码 ↔ Table 4-4 / 附录 A 的对应 | Table 4-4、附录 A |
| `notebooks/01_replication_walkthrough.ipynb` | 带说明的分步运行笔记本 | — |
| `outputs/` | 脚本生成的结果表与日志，可随时删除重跑 | — |

## 二、上传到 GitHub 的步骤

1. 在本地解压本材料包，进入 `Research-data-Teesside-PhD-thesis/` 目录。
2. 先在本地跑一遍，确认环境正常：
   ```bash
   python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   python code/run_all.py
   ```
3. 初始化并推送：
   ```bash
   git init
   git add .
   git commit -m "Replication materials: code, Mplus syntax, synthetic data (thesis v5.2)"
   git branch -M main
   git remote add origin https://github.com/xianzifeng/Research-data-Teesside-PhD-thesis.git
   git push -u origin main
   ```
   （若仓库已存在且非空，先 `git pull origin main --allow-unrelated-histories`。）
4. 在 GitHub 上 **Releases → Draft a new release**，tag 填 `v1.0.0`，标题 `Thesis submission version`。
5. 获取 DOI：登录 <https://zenodo.org>（用 GitHub 账号），在 *GitHub* 页面把本仓库的开关打开，然后回到 GitHub 发布 release；几分钟后 Zenodo 会生成 DOI。
6. DOI 已生成并回填：`10.5281/zenodo.23448040`（release 1.0.1，2026-10-18）。README 徽章、`CITATION.cff` 与论文 4.8.4、4.9.4、附录 C、H.4、Table 7-2 均已引用该 DOI。以后若发新 release，Zenodo 会生成新的版本 DOI；论文引用保持 1.0.1 的版本 DOI 即可。

## 三、答辩前还需作者补充的内容

- `CITATION.cff` 中的 ORCID（可选）。
- 若实际使用的 Mplus 输入与此处模板有差异（例如控制变量进入方式、二阶结构、残差相关），请以**实际运行过的 .inp** 替换，并把 `.out` 文件一并提交到 `mplus/outputs/`——这是审稿人评论 C776 明确要求的。
- 若希望保留真实分析的结果表，可把 Mplus/Python 在真实数据上的输出（不含个体数据）放入 `outputs/real_results/`。
- 模拟数据的 IV 第二方程（RR 被工具变量化）第一阶段 F 较弱，这是模拟数据的性质，与论文结论无关；真实数据运行时以实际输出为准。

## 四、模拟数据的生成原理（供答辩时解释）

模拟数据**不含任何真实观测**。生成器只使用论文中已公开的参数：Table 5-1 的均值、标准差与相关系数；Table 4-4 的标准化载荷；Table 5-3 / 5-13 的标准化路径（PR1→FA 0.412、PR2→FA −0.198、PR1×PR2 −0.124、FA→RR 0.399、RR→DC 0.452、FA→DC 0.267、PR→DC 0.134）；§4.9.3 的员工层 ICC（0.42 / 0.38 / 0.29）；Table 5-6/5-7 的 ET、DA 分组比例。先按结构方程生成标准化潜变量，再按载荷生成 1–5 点量表题项并取整。`data/synthetic/generation_log.txt` 列出每个目标值与实际达到值的对照。
