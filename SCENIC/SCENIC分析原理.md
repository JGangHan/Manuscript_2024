**SCENIC (Single-Cell rEgulatory Network Inference and Clustering) 目前其实只有两个在正常维护，分别是 pySCENIC 和 SCENIC+。另外的SCENIC（R）和VSN-pipeline（就一个流程指南，看了看没什么有用的东西）目前已经不推荐使用了**
![image](https://github.com/JGangHan/sc_data_analysis/assets/75400599/58c998eb-20b0-46ae-8272-0a37ddc6a758)
# pySCENIC
**是基于 GENIE3，cisTarget和AUCell三个R包构建的**
## Tips
1. 虽然是做转录因子分析，但其实前边也有细胞质控步骤，跟Seurat质控步骤相似。**这个软件里边的折线图是比较好看的**
2. 细胞质控后也是PCA，高变基因鉴定，clustering
3. 再之后是转录因子分析
## 分析原理
1. **CO-expression**. 用**GENIE3**软件将转录因子和其他基因进行共表达分析，鉴别到众多TF-Genes共表达模块（module/regulon），一个转录因子对应一个模块，不同模块可以有共同的转录因子
2. **Motif discovery**。用**cisTarget**数据库，这个数据库中包含人/小鼠转录因子motif序列，检测每个TF-Genes共表达模块中每个基因上游是否包含该转录因子motif序列，过滤掉不包含该转录因子motif序列的基因
3. **Cell scoring**. 计算AUCell algorithm用来表示每个模块在单个细胞中的 regulon activit，从而构建一个 AUC score x cells 矩阵（类似genes x cells矩阵），从而可以用于数据降维
![image](https://github.com/JGangHan/sc_data_analysis/assets/75400599/62da89c5-e2f4-471f-beb8-5c4970e7d89c)










