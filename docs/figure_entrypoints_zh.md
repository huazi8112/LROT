# 论文图件运行入口

[English](figure_entrypoints.md)

在仓库根目录、记录环境中执行。本表列出支持的主要入口，不枚举全部历史辅助脚本。将 `<out>` 替换为仓库外尚不存在的新目录，每次调用使用不同目录；将 `<S>` 替换为该行列出的一个补图编号。字体必须是哈希与清单匹配的本地文件。原始输入及外部依赖需另行准备并匹配身份，命令不自动下载。

先核验发布文件；需要论文正式图片时，再使用原样导出命令。以下命令不运行实验：

```text
python tools/verify_manifest.py
python figures/render_current.py --mode asset --output-dir "<正文图片新目录>"
python figures/render_supplement.py --mode asset --output-dir "<补充图片新目录>"
```

下表输入列仅列主要输入，完整执行绑定见[正文清单](../figures/current_release/figure_manifest.json)和[补图清单](../figures/supplement/figure_manifest.json)；正式图片身份见[31图索引](../figures/current_release/FINAL_PAPER_IMAGE_INDEX.csv)。部分源码文件名沿用旧图号，应以本表的现行图号为准。

| 图号／任务 | 命令 | 主要输入 | 验证范围 |
|---|---|---|---|
| F1 | `python figures/render_current.py --figure F1 --mode records --output-dir "<out>"` | [F1.png](../figures/current_release/assets/F1.png) | 仅重绘 LR 先验标签区域，保留其他像素；示意图操作，不运行实验。 |
| F2 | `python figures/render_current.py --figure F2 --mode records --output-dir "<out>"` | [80步保存曲线](../figures/current_release/inputs/F2_four_gamma_histories_current_code.npz); [F2.png](../figures/current_release/assets/F2.png) | 由保存曲线重绘 F 面板；保留 A–E 和123值主轨迹，不认证历史恢复或收敛。 |
| F3 | `python figures/render_current.py --figure F3 --mode records --output-dir "<out>"` | [F3.json](../figures/current_release/inputs/F3.json) | 绘制保存指标及配对统计；STAligner 评分/训练、30组求解见本表后续行。 |
| F4 | `python figures/render_current.py --figure F4 --mode records --output-dir "<out>"` | [DLPFC NPZ](../lrot_output/arrays/lrot_data_dlpfc.npz); [保存汇总](../lrot_output/dlpfc_results/lrot_dlpfc_results.txt) | 仅由保存数组和结果绘图，不重新执行原始 DLPFC 实验。 |
| F5 | `python figures/render_current.py --figure F5 --mode records --output-dir "<out>"` | [逐 spot 数组](../figures/current_release/inputs/F5_current_protocol_plot_only.npz); [输入绑定](../figures/current_release/inputs/F5_input_binding.json) | 绘制保存的逐 spot 结果；重绘并非正式图字节完全一致，正式图片应使用 asset 导出。 |
| F6 | `python figures/render_current.py --figure F6 --mode records --output-dir "<out>"` | [鼠脑数组](../lrot_output/arrays/lrot_data_real_visium.npz) | 保存结果绘图；实测表达配合生成坐标和扰动，区域准确率不适用。 |
| F7 | `python figures/render_current.py --figure F7 --mode records --output-dir "<out>"` | [历史计划绘图输入](../figures/current_release/inputs/F7_historical_plan_plot_only.npz); [绑定](../figures/current_release/inputs/F7_source_binding.json) | 由保留计划及后续匹配掩膜汇总绘图；历史同次运行来源未绑定，新固定输入运行另列。 |
| F8 | `python figures/render_current.py --figure F8 --mode records --output-dir "<out>"` | [F8.json](../figures/current_release/inputs/F8.json); [三个面板](../figures/current_release/components/) | 合成三个保存面板，不重新执行三维对齐或多种子实验。 |
| S2, S4, S5, S7, S8, S9, S15, S17, S19 | `python figures/render_supplement.py --figure "<S>" --mode records --output-dir "<out>"` | [补图清单](../figures/supplement/figure_manifest.json); [保存输入](../figures/supplement/inputs/); [S19 计划与绑定](../figures/supplement/s19_posthoc/) | 在批准画布上重绘保存结果；S4A 使用报告精度的汇总，S19B 从保存计划重建；保留面板不等于独立重算。 |
| S3, S6, S12, S20, S21 | `python figures/render_supplement.py --figure "<S>" --mode patch --font-regular "<arial.ttf>" --font-bold "<arialbd.ttf>" --output-dir "<out>"` | [局部修订记录](../figures/supplement/inputs/); [字体绑定哈希](../figures/supplement/figure_manifest.json) | 使用绑定字体进行局部结果/标签修订，不运行新实验，不代表整图由原始输入重建。 |
| S1, S10, S11, S13, S14, S16, S18, S22, S23 | `python figures/render_supplement.py --figure "<S>" --mode asset --output-dir "<out>"` | [正式补图](../figures/supplement/assets/) | 本主入口提供正式图片原样导出；S23 历史点位移修补入口已停用。 |
| F2：新诊断曲线 | `python figures/reconstruct_f2_curves_current.py --output-dir "<out>"` | [模拟输入归档](../lrot_output/arrays/lrot_data_synthetic.npz); [固定设计](../experiments/synthetic_design_manifest.json) | 执行四条固定80步诊断曲线；123值主轨迹仍为保留记录，不证明收敛或恢复历史运行。 |
| F3：STAligner 评分 | `python experiments/reproduce_f3_staligner.py --mode verify --output-dir "<out>"` | [有序输入和两次保存嵌入](../provenance/f3_staligner/) | 仅评分保存嵌入，不训练；固定案例一致率 .930、ARI .853503，不是旧 .934 的恢复。 |
| F3：STAligner 训练 | `python experiments/reproduce_f3_staligner.py --mode train --output-dir "<out>"` | [固定案例与环境绑定](../provenance/f3_staligner/README.md); [离线源码准备](external_resources.md); [依赖绑定](../provenance/f3_staligner/environment.json) | 一次完整、确定性的1000步固定案例训练；需匹配外部源码与依赖。不自动下载，不搜索种子，不认证全文实验。 |
| F3C / S23：30组配对 | `python experiments/supp/synthetic_seed_robustness.py --verify-only`<br>`python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<out>"` | [配对注册表和汇总](../lrot_output/s23_support/); [绑定数值源码](../paper_protocol/source_binding.json) | 第一条只读核对保存记录；第二条在新目录、不加 --resume 时执行30组、60次新求解。历史记录与当前结果分别保留身份。 |
| F5：Br5292 原始输入案例 | `python experiments/validate_f5_protocol.py --repo-root . --data-dir "<dlpfc-dir>" --output-dir "<out>"` | [冻结输入绑定](../figures/current_release/inputs/F5_input_binding.json); [外部 DLPFC 文件](../data/README.md) | Br5292 固定案例四次求解及独立逐 spot 评分；不恢复历史运行绑定，不认证第二供体。 |
| F6：鼠脑原始输入案例 | `python experiments/run_real_visium.py --solver revalidated --input-h5 "<mouse.h5>" --output-dir "<out>"` | [外部鼠脑 H5及校验值](../data/README.md); [输入顺序证据](../provenance/mouse_perturbation/) | 按记录精度和绑定求解器重算固定表达/扰动案例；坐标为生成坐标，无区域准确率。 |
| F7：新固定输入案例 | [四步命令：准备→求解→绘图输入→绘图](../provenance/f7_prospective/README.md) | [原始文件和顺序绑定](../provenance/f7_prospective/reference_input_identity.json); [外部乳腺癌文件](../data/README.md) | 核对原始输入后进行两次新求解，评分/绘图使用同次新计划；正式历史图另行保留。准备目录须已存在且为空，后续阶段使用新目标。 |

图像导出、保存结果重绘、保存嵌入评分及新数值执行的验证范围不同。绘图一致不能认证全部实验条件。环境准备见[运行指南](../RUN_REPRO_zh.md)，F7 多阶段目录要求见[案例说明](../provenance/f7_prospective/README.md)，证据边界见[验证范围](release_readiness.md)及[数值验收记录](../provenance/reproduction_acceptance/summary.json)。本次只更新导航文档，不执行实验或替换论文图片。


## Targeted checks

| Paper item | Entry | Scope |
|---|---|---|
| Tables S9/S22 | `experiments/prepare_dlpfc_staligner.py` → `experiments/train_dlpfc_staligner.py` | Two explicitly input-bound fresh STAligner cases; preparation and training use separate processes. |
| External comparisons and GO | `python experiments/verify_comparators.py --output-dir "<new-out>"` | Independent scoring of 106 retained plans, four historical/fresh DLPFC embeddings and 14,584 GO term statistics. No optimiser or training. |
| Figure S3 statistics / Note S16 | `python experiments/go_enrichment.py --obo "<go-basic.obo>" --gaf "<mgi.gaf.gz>" --output-dir "<new-out>"` | Regenerates the complete term table from exact resource checksums; approved image bytes remain unchanged. |

Commands, inputs and evidence levels: [targeted verification guide](comparator_verification.md).
