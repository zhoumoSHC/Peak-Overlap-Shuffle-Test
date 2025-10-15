# Peak-Overlap-Shuffle-Test
A randomization-based statistical framework to evaluate whether the observed overlap between two sets of genomic regions (peaks) is significantly higher than expected by chance.
这个分析流程用于评估两个基因组区间（BED 文件）之间的重叠是否显著高于随机期望。

**一、运行随机化检验**

在终端运行：
```
  bash peak_overlap_shuffle_plot.sh \
  -a A.bed -b B.bed -g hg38.chrom.sizes \
  -w 150 -n 1000 -o overlap_results.txt \
  --seed 20251015
```
参数说明：

`-a`, `-b`：输入的两个 BED 文件

`-g`：基因组大小文件（如 hg38.chrom.sizes）

`-w`：重叠窗口大小（bp）

`-n`：随机化次数

`-o`：输出随机结果文件

运行结束后会输出实际重叠数、均值、Z-score、p-value。

**二、绘制随机分布图**

在 R 中运行“shuffle作图.R”，注意修改overlap的数值和输出的文件名。
