# LROT

本目录可独立作为代码仓库根目录。固定GO/MGI输入已纳入`data/go/`，S1工作簿
已纳入`supplementary_files/`，不依赖外层交付目录或另一个未发布的资源包。
见[仓库整理与检查说明](docs/repository_upload.md)。

用于空间转录组切片对齐的配体–受体引导最优传输方法。

LROT 将表达相似性、人工指定的配体–受体（LR）兼容性先验和切片内部空间
结构结合起来，求得两张切片之间的软传输计划。包内保留实现、固定执行协议、
所选结果、绘图输入和核验记录。LROT 求解器本身不需要神经网络训练。

[English](README.md) · [运行说明](RUN_REPRO_zh.md) ·
[逐图入口](docs/figure_entrypoints_zh.md) · [核验范围](docs/release_readiness.md)

## 目录与输入

| 位置（相对代码目录） | 内容 |
|---|---|
| `paper_protocol/` | 源码哈希校验、数值实现绑定及固定求解包装器 |
| `provenance/numerical_validation_snapshot/code/` | 固定协议所用的保留计算源码 |
| `experiments/` | 分实验的数据准备、求解和评价入口 |
| `figures/` | 所选原图、保存绘图输入、清单和绘图入口 |
| `lrot_output/` | 所选汇总结果、保存数组和种子配对记录 |
| `data/`、`cancer_data/`、`dlpfc_data/` | 数据来源、文件校验、准备步骤与许可说明 |
| `docs/`、`tools/` | 环境、方法运行说明和完整性检查工具 |
| `data/go/` | 固定 GO 本体、注释、来源和校验记录 |
| `supplementary_files/` | Supplementary File S1 工作簿 |

原始组织矩阵和 CellChatDB 二进制另行取得。压缩包不包含论文 DOCX/PDF、
虚拟环境、凭据、上传工具、临时日志或私有备份。第三方数据和源码的来源链接
用于识别所需资源，不是本代码包的托管地址。

## 完整方法流程

### 1. 固定输入身份和预处理

准备切片 A/B 的有序 spot×gene 表达数组、相应坐标和共同的有序基因列表。
样本身份、条码、坐标单位、变换和预处理由具体实验的输入绑定确定；层标签仅
用于需要注释的评价。原始计数、转换后的表达和保存绘图数组不是同一种输入。

严格使用相应加载器和已记录配置。保持抽样索引原顺序、数据精度、随机种子、
配对方向及基因顺序，不能为匹配某个结果重新排序或更换归一化。
[数据说明](data/README.md)和[输入清单](data/input_files.json)给出来源及 SHA-256。
DLPFC 以切片编号识别；乳腺癌使用 Block A 的 Section 1/2。小鼠扰动案例是一份
实测表达数据配合生成坐标和施加变换，不是两张独立实测连续切片。
模拟数据遵循固定设计和已记录种子；registered 在这里表示记录在案，不表示预注册。

### 2. 选择物种对应的 LR 先验

使用实现中的小鼠或人类配对及原权重。权重是人工兼容性先验，不是实测作用
强度，也不是数据库置信分数。按原实现跳过缺失基因；基因符号处理和表达掩膜
以具体入口为准。不能替换物种表，也不能仅凭大小写转换认定同源关系。

逐对资源含义与限制见[LR 说明](docs/lr_resources.md)和
[证据记录](data/lr_resource_evidence.json)。
[S1 工作簿](supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx)保留鼠/人表。

### 3. 构建表达、LR 和几何三个通道

- **表达通道：** A/B 表达谱之间的余弦距离，按代码中的最大值与数值保护归一化。
  默认表达掩膜排除 LR 基因；当非 LR 特征为空或相关基因信息缺失时，按源码
  保留的回退规则使用全基因，不自行更换掩膜。
- **LR 通道：** 对可用配对计算
  `S_raw[i,j] = sum(w[L,R] * expression_A[i,L] * expression_B[j,R])`，
  再作 min–max 归一化。常量矩阵置零。非恒定 S 的代价为 `-log(S + 1e-4)`，
  再按实现归一化；常量 LR 通道代价为零。
- **几何通道：** 计算各切片内部的平方欧氏距离并分别按最大值归一化。
  Gromov–Wasserstein 代价包含 `-2 * D_A @ P @ D_B.T`，每次外层更新都随当前
  P 重新计算并缩放，不能把它视为固定的切片间坐标距离。

精确实现见[绑定的核心源码](provenance/numerical_validation_snapshot/code/lrot_core.py)。
代码根目录的 `lrot_core.py` 用于明确标记的旧协议，不能替代当前绑定入口。

### 4. 按固定数值协议求解

`paper_protocol.binding.bind` 校验固定源码并显式绑定经数值验收的平衡求解器；
`solve_locked` 在运行前核对预期输入身份。入口和哈希见
[binding.py](paper_protocol/binding.py)及[源码绑定清单](paper_protocol/source_binding.json)。

采用均匀边际，初始计划为外积 `P = a[:,None] * b[None,:]`。每次更新使用：

```text
C_total = alpha * C_expression + gamma * C_LR
          + (1 - alpha) * beta * C_GW(P)
effective_epsilon = reg_ot + beta * reg_gw
```

内层先进行最多1,000步缩放，必要时再进入缩减对偶优化和有界可行性修正。
是否接受由残差和证书检查共同决定，不能只看优化器的 success。
包装器中的 `inner_max_iter=10000` 不表示初始缩放实际执行10,000步。

固定包装器采用边际绝对容差 `1e-8`、相对容差 `1e-5`，外层相对计划变化容差
`1e-5`；最早第7次外层更新后接受。只允许包装器原定的外层上限重试，失败
结果不得进入图表。具体实验的参数、上限及尝试记录仍是该实验的依据。

保存的 `loss` 是代码中的诊断代价汇总，并不是完整的熵正则更新目标；不能
仅凭曲线宣称全局最优或收敛。不能换近似算法、种子、精度或归一化来追求旧值。

### 5. 重心映射与连续堆叠

A→B 映射对 P 按行归一化，计算
`mapped_A = (P / row_sums(P)) @ coordinates_B`，并使用相应入口的数值保护。
反方向映射需要列归一化及转置计划；两者的 spot 身份和方向不同。

解释计划前应核对 shape、有限性、非负性、边际残差和零质量行/列。
连续堆叠按记录的配对、变换和顺序执行；一对切片通过不代表整个三维实验通过。

### 6. 按原指标和方向评价

| 指标 | 定义及边界 |
|---|---|
| soft 标签转移准确率 / ARI | 行归一化计划将 B 的已知层标签转移到 A，再与 A 标签比较；依赖注释，不是对应真值 |
| hard 传输匹配 | 每行取最大质量的 B spot 后比较标签，与 soft 转移不同 |
| 对齐坐标最近邻准确率 | A 重心映射后，在原始 B 坐标中寻找最近邻并比较标签，与最大质量匹配不同 |
| 行熵 | 行归一化计划的平均熵；较低表示更集中，不直接证明更准确 |
| 同标签传输质量 | 同标签连接所占质量比例，与逐 spot 准确率不同 |
| 分实验评分与汇总 | 保持原评分、掩膜、分母、排序及聚合；不能从均值反推出逐 spot 数组 |

不得把不同指标混为统一排名，也不得将注释派生分数写成标定概率。
STAligner 的嵌入评分属于另一条流程，其完整训练和依赖分别固定。
计时需说明准备、求解和评分的边界、机器及重复次数；单次本地观察不能支撑普遍速度排名。

### 7. 保存、核验与绘图

记录原文件和有序数组身份、参数、种子、源码哈希、环境、完整尝试和诊断、P、
指标及日志。新输出使用包外未存在的目录；保留历史结果，并明确新运行与历史记录。

区分原图资产导出、由保存结果绘图和重新数值计算。
[逐图入口表](docs/figure_entrypoints_zh.md)及图件清单注明每张图支持的层级。
PNG 相同不等于实验端到端复现。

## 环境与快速开始

记录环境为 Windows 11 / AMD64、Python 3.10.11。准备方式见
[已测试环境说明](docs/clean_environment.md)和[分任务依赖](docs/environments/)。
[requirements.txt](requirements.txt)是较宽的研究依赖清单，不能据此认定所有外部
方法已就绪。入口不会自动安装依赖。

在代码目录中执行以下命令。输出目录必须尚不存在；下例仅导出原图，不求解：

```powershell
python tools/verify_manifest.py --strict
python tools/verify_inputs.py --list
python figures/render_current.py --list
python figures/render_supplement.py --list
python figures/render_current.py --mode asset --output-dir "..\..\lrot_main_assets"
python figures/render_supplement.py --mode asset --output-dir "..\..\lrot_supp_assets"
```

由保存记录重画 F2 支持的面板，而不重新求解：

```powershell
python figures/render_current.py --figure F2 --mode records --output-dir "..\..\lrot_F2_records"
```

该入口保留 A–E 和单独的123值主曲线，利用保存诊断记录重画 F 面板。
不同图件的支持范围不同；精确像素一致还依赖绘图环境与字体。

## 明确选择的数值运行入口

下列命令仅供有意开展数值运行时使用，需要对应数据和环境；本次打包没有执行：

```text
python experiments/supp/synthetic_seed_robustness.py --verify-only
python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<new-30-pair-output>"
python experiments/validate_f5_protocol.py --repo-root . --data-dir "<DLPFC-directory>" --output-dir "<new-F5-output>"
python experiments/run_real_visium.py --solver revalidated --input-h5 "<mouse-H5>" --output-dir "<new-mouse-output>"
```

30组模拟运行在不复用检查点时执行60次新求解。F3/STAligner、F7有独立的准备、
执行和评价步骤，见[运行说明](RUN_REPRO_zh.md)和
[F7流程](provenance/f7_prospective/README.md)。不要直接启动旧 driver 的 main
来替代绑定协议。

## 附件、资源与复现边界

S2是[一份包含30组种子配对的CSV](lrot_output/s23_support/Supplementary_File_S2_synthetic_seed_pairs.csv)，
不是30份CSV。其他实验的种子见各自配置和证据。
[固定GO资源](data/go/README.md)附带来源、许可和校验记录；组织数据及
CellChatDB二进制另行准备。MIT代码许可不替代第三方资源条款，见
[外部资源说明](docs/external_resources.md)。

文件完整性、原图一致性、保存结果绘图、保存数组独立评分和从原始输入重新
求解是不同层级。具体证据与限制见[核验范围](docs/release_readiness.md)。
并非所有论文实验均已端到端复现；技术重复不能代替生物学重复。
保留历史案例标识用于追溯。

## 许可与快照标识

源码许可见[LICENSE](LICENSE)。用[LOCAL_VERSION.json](LOCAL_VERSION.json)、
[RELEASE_FILES.txt](RELEASE_FILES.txt)及校验清单识别本地快照，外层压缩包另附 SHA-256。
[CITATION.cff](CITATION.cff)保留既有集体贡献者条目；本包不代填个人作者、
公开托管地址、正式软件发布版本或归档 DOI。
