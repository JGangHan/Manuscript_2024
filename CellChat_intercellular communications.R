library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

### 1. extrace gene expression information
# removal of unknow cells
data.input = levels(ob.merge.rmdb$pie_celltype)
data.input = data.input[-23]
data.input
data.input = subset(ob.merge.rmdb, pie_celltype %in% data.input)
# extract normalized data
data.input <- GetAssayData(data.input, assay = "RNA", slot = "data")
head(data.input)
row.names(data.input)

# transfer sheep gene symbol to human homologous gene 
name = row.names(data.input) %>% as.data.frame()
colnames(name) = "id1"
# 原始 id 转为 ensembl id
name$id2 = id_map_seurat[name$id1, 1]
# ensemble id 转为 merged_Symbol_deoverlap
name$id3 = id_map_v2[name$id2, 7]
# 检验基因id顺序是否一致
identical(row.names(data.input), name$id1)
# 重新赋值 data.input 行名
row.names(data.input) = name$id3
row.names(data.input)


### 2. cell meta information
data.input = levels(ob.merge.rmdb$pie_celltype)
data.input = data.input[-23]
data.input
data.input = subset(ob.merge.rmdb, pie_celltype %in% data.input)
meta = data.input@meta.data


### 3. create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype_final")
data.input = NULL
# select huamn ligand-receptor interaction database
CellChatDB <- CellChatDB.human # CellChatDB.mouse for mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use


### 4. interaction pathways enrichment analysis
# select pathway-related genes
cellchat <- subsetData(cellchat, features = NULL) 
# identification of overexpressed genes
cellchat = identifyOverExpressedGenes(
  cellchat, thresh.pc = 0, thresh.fc = 0.5, thresh.p = 0.01,
  data.use = NULL,    # 可以用自主构建的数据，一般用不到
  group.by = NULL,    # cell group information, one column names of meta dataset
  idents.use = NULL,  # 用来分析的细胞子集
  group.dataset = NULL, # 当涉及多个合并 cell chat 对象时设置该参数
  pos.dataset = NULL, # 多个合并 cell chat 对象中分析目标 cell chat 对象名称
  features.name = "features" # 计算结果储存名称
)
# OverExpressedInteractions
cellchat <- identifyOverExpressedInteractions(cellchat, return.object = T) # 无其他参数
# 去掉不需要的 levels
levels(cellchat@idents)
cellchat@idents = droplevels(cellchat@idents)
levels(cellchat@idents)
# 计算 LR pairs 水平细胞互作概率
cellchat <- computeCommunProb(
  cellchat, type = "triMean", # "triMean", "truncatedMean", "thresholdedMean", "median"
  # 如果设置 type = "truncatedMean"，需要指定 trim 值
  trim = 0.1,  
  LR.use = NULL,        # 指定 LR paris
  distance.use = TRUE,  # 是否计算互作概率的空间距离
  interaction.length = 200, # 配体的最长互作距离
  scale.distance = 0.01, # 比较两组不同的单细胞数据时设置该参数
  k.min = 10,  # the minimum number of interacting cell pairs required for defining adjacent cell groups
  nboot = 100 # p 阈值
)
# LR pairs
# 过滤 低表达细胞
cellchat <- filterCommunication(cellchat, min.cells = 10) 
# 计算通路水平细胞互作
cellchat <- computeCommunProbPathway(cellchat)


### 5. final 297 genes, 143 LR pairs, 21 pathways
cellchat@var.features$features %>% length()
cellchat@net$prob %>% dim() 
cellchat@netP$prob  %>% dim() 


### 6. merge interaction networkss
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, './cell_communication/cellchat.rds')


### 7. save data
# 所有 LR pairs 互作
df.net_LR <- subsetCommunication(cellchat, slot.name = "net")
# 共 149 LR pair
table(df.net_LR$interaction_name) %>% length()
# 所有 通路水平 互作
df.net_Path <- subsetCommunication(cellchat, slot.name = "netP")
# 149 LR pair 分配在 22 个通路中
table(df.net_Path$pathway_name) %>% length()
write.csv(df.net_Path, './1.csv')
write.csv(df.net_LR, './1.csv')










