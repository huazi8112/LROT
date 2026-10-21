# 运行与复现

先按[图件入口表](docs/figure_entrypoints_zh.md)选择图号和操作，再按下文准备环境、输入并确认验证范围。

在仓库根目录运行。使用 Python 3.10.11。新安装按[已验证任务依赖及锁定清单](docs/clean_environment.md)执行；`requirements.txt` 保留较广的研究环境记录，补齐了此前漏列的两项 PyG CPU 扩展，已有数值依赖版本不变。STAligner 完整训练另需匹配 `provenance/f3_staligner/environment.json` 的外部源码和依赖；以下工具不自动安装或下载。

## 1. 准备输入

[data/README.md](data/README.md) 统一列出数据编号、获取地址、文件路径及校验范围。原始组织输入及CellChatDB二进制不随代码包提供；固定GO/MGI文件已在`data/go/`。外部组织数据根目录可含`cancer_data/`、`dlpfc_data/`：

```text
python tools/verify_repository_package.py
python tools/verify_inputs.py --dataset go --data-root .
python tools/verify_inputs.py --list
python tools/verify_inputs.py --dataset breast --dataset dlpfc --data-root "<数据根目录>"
python tools/verify_inputs.py --dataset cellchat --resource-dir "<已授权资源目录>"
python tools/verify_manifest.py
```

只选择运行需要的数据集；缺失或哈希不匹配即报错。新增盘点哈希标识当前本地文件，不补造历史获取/冻结记录。CellChatDB 再分发依据仍待确认，需另行提供有权使用的匹配文件。

## 2. 导出或重绘图件

输出目录必须位于仓库之外且尚不存在，请替换占位符：

```text
python figures/render_current.py --list
python figures/render_supplement.py --list
python figures/render_current.py --mode asset --output-dir "<正文资产新目录>"
python figures/render_supplement.py --mode asset --output-dir "<补充资产新目录>"
python figures/render_current.py --figure F3 --mode records --output-dir "<F3重绘新目录>"
python figures/render_supplement.py --figure S4 --mode records --output-dir "<S4重绘新目录>"
```

`asset` 原样导出8张正文图或23张补图，不求解；`records` 用保存输入重绘已支持图件，不等于重新执行实验。正文默认采用 `records`。入口核验文件哈希；补图模式见[补图说明](figures/supplement/README.md)，选图位置见[图件索引](lrot_output/FIGURE_INDEX.md)。

S4A 使用保存的四位小数汇总，B/C 像素保留；S19B 使用保存计划及有序坐标/标签，直接以当前批准图保留 B 面板重绘区域以外的像素，不再保留旧模板。字体和环境可能影响绘图字节；重绘失败不能自动改写论文批准图。

## 3. 数值案例复现

修订数值执行使用 `paper_protocol/`；显式 legacy 模式仅用于历史协议。输入、精度、种子及求解器绑定均保持固定。

```text
python experiments/supp/synthetic_seed_robustness.py --verify-only
python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<30组复现新目录>"
python experiments/reproduce_f3_staligner.py --mode verify --output-dir "<F3评分新目录>"
python experiments/reproduce_f3_staligner.py --mode train --output-dir "<F3训练新目录>"
python experiments/run_real_visium.py --solver revalidated --input-h5 "<鼠脑H5>" --output-dir "<鼠脑新目录>"
python experiments/validate_f5_protocol.py --repo-root . --data-dir "<DLPFC目录>" --output-dir "<F5新目录>"
```

30组 `--run` 执行60次新求解，不复用检查点。F3 的 `verify` 仅评分保存嵌入；`train` 完整执行确定性训练，外部源码或依赖绑定变化时报错。F5 对 Br5292 案例执行四次求解。鼠脑采用实测表达、生成网格坐标及人工扰动，无区域准确率；既定 float16 转 float32 输入协议不变。

F2 曲线命令：`python figures/reconstruct_f2_curves_current.py --output-dir "<F2新目录>"`。计算四条80步诊断曲线；另存的123值主轨迹仍是保留记录，不宣称新运行恢复。F7 输入准备、求解及评分绑定见[案例说明](provenance/f7_prospective/README.md)。

## 4. 检查范围与发布

[docs/release_readiness.md](docs/release_readiness.md) 集中说明现行证据和限制；新数值验收保存在 `provenance/reproduction_acceptance/summary.json`。文件完整性、图像一致和技术重复均不等于全部实验端到端认证或生物学重复。

`RELEASE_FILES.txt` 定义分发范围。用 `python tools/build_release.py --output "<新发布包目录>"` 构建，再在该副本运行 `tools/verify_manifest.py --strict`。Supplementary File S1已纳入`supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx`；所选30行Supplementary File S2 CSV也已随包。附件编号与 Figure S1/S2 图号不同。新运行输出及原始资源不自动纳入发布。

[English](RUN_REPRO.md)

## Fixed external-source preparation

See [external_resources.md](docs/external_resources.md) for pinned public CellChatDB acquisition/exact PKL rebuilding and offline reconstruction of the six fixed STAligner training files. The original training source/version/backend gates remain enforced. This source recipe does not certify a new clean dependency installation or recover historical .934.


定向核验：两对 DLPFC 的新训练结果为准确率 0.586/0.525，历史嵌入重新评分仍为 0.580/0.518，两套记录分别保留；106 份外部比较计划已核对输入依据并独立评分；固定 GO/MGI 原始文件逐条重建全部 14,584 个检验结果，并与官方当前下载逐字节匹配。入口、命令和边界见[定向核验说明](docs/comparator_verification.md)。这些核验不等于全部论文实验均已端到端复现。
