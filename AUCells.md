## [tutorail](https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#before-starting)
## [github code](https://github.com/aertslab/AUCell/tree/master/R)
**tips：**
由于AUCells在分析过程中会将原始gene * cell表达矩阵转换为 ranking * cell排序矩阵，所以在原始的表达矩阵中的数值无论是raw count，TMP，UMI count或者经过归一化方法矫正的数值，都不会影响后续的分析结果，这种分析方法理论上避免了library size对基因集活性分析的影响

## 分析流程
1. **count matrix**转为sparse count matrix；多个**gene sets**，这些基因集可以包含数量不等的基因

2. 根据count matrix中每个基因的表达量，利用colRanks()函数对每一个细胞中所有的基因进行排序，输出结果是**一个细胞中所有基因的表达量排名顺序**，也就是将表达量矩阵转化为了排名矩阵 **ranking matrix**，用于后续分析

3. 计算恢复曲线（revovery curve）的“Area Under the Curve” (AUC)。
   * 首先是怎么理解recovery curve：恢复曲线是一种用于评估基因集在基因排名中的富集情况的图形工具。具体来说，它是一个绘图，其横坐标表示从排名最高的基因开始累积考虑的基因数量（从0个基因开始逐渐增加，到10，50，100），纵坐标表示这些基因中属于特定基因集的基因的比例（也是从0%开始逐渐增加，但绝对到不了100%）。换句话说，这条曲线展示了在考虑越来越多的基因时，特定基因集的基因如何逐渐被“恢复”（即被识别出来）。
   * 需要注意基因集被恢复的比例绝对到不了100%，因为为客观评价基因集活性，会在恢复曲线中设定阈值，例如仅观察ranking matrix中排名前5%的基因对基因集的恢复情况。这个阈值是可以修改的
   * 基因集被恢复的比例越大，表明这个基因集中越多的基因整体排名靠前，基因表达水平越高，基因集活性越高
   * 当计算多个基因集的AUC score, 就会生成一个**AUC scores*cells矩阵**
   * **需要注意的是最初步得到的AUC socre是绝对值，还可以选择是否进行额外的归一会处理，进而使其在不同实验或不同条件下具有可比性**
![image](https://github.com/JGangHan/sc_data_analysis/assets/75400599/877a56eb-dbf4-4383-8885-5d1041444990)

4. AUC score 分布特征：因为多种细胞类型和多个基因集的存在，可以探究指定基因集AUC socre在所有细胞中的密度分布特征。
   * 双峰分布（bimodal）：如果一个或少数几个基因集在特定细胞类型中非常活跃，AUC socre高，此时出现一个频率顶点；而其他基因集在这些细胞中表达低，出现第二个频率顶点，从而呈现为双峰形态。此时选择双峰之间的拐点作为阈值；
   * 长尾分布（distribution with a long tail）：这种分布特征的出现表明某些基因集非常活跃，但不是在全部细胞中或绝大多数细胞中，从而出现长尾，AUC score越高的基因集，出现频率越低。此时选择尾部离群值作为阈值；
   * 正态分布（bimodal）：通常为非特异表达基因集分布特征，比如随机基因集或管家基因集，理论上这些基因集没有进一步分析意义，但需看具体情况而定
   * ![image](https://github.com/JGangHan/Software-list/assets/75400599/c6c6ba3b-1f10-4841-b32c-8e8b1b2ef4df)

5. AUC score阈值选择：
   * 这一步的目的是如何设置一个合理的AUC阈值，从而筛选AUC score大于该阈值的基因集，聚焦具有更重要生物学意义的细胞集
   * 可以通过 **AUCell_exploreThresholds()** 函数自动设定或手动设定

6. 后续分析
   * 观察指定基因集的AUC score在细胞中出现频率<img src="https://github.com/JGangHan/sc_data_analysis/assets/75400599/da75219b-f050-4bfb-9071-81e090c3fea4" width="500" height="300">
   * 利用阈值将AUC score*cell matrix转为二进制矩阵<img src="https://github.com/JGangHan/sc_data_analysis/assets/75400599/efd0e2d3-7267-457b-9cbd-187ff52e9ef2" width="500" height="300">
   * 用AUC score对细胞tSNE/UMAP降维进行可视化<img src="https://github.com/JGangHan/sc_data_analysis/assets/75400599/8ef0c944-0cef-4e5b-9b22-b76822892615" width="500" height="300">




