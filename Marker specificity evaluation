########################## example cell type-adipocyte ############################

library(Seurat)
library(tidyverse)
library(dplyr)


### 1. DEGs analysis
marker_adipocyte = FindMarkers(
  ob.merge.rmdb, logfc.threshold = 0.5, min.pct = 0, only.pos = T,
  ident.1 = c("Adipocyte"))
data = marker_adipocyte


### 2. gene annotation
# gene
data$Gene = row.names(data)
# ensembl
data$Ensembl=id_map_seurat[data$Gene,"Ensembl"]
# gene_merged
data$Gene_merged = id_map[data$Ensembl, "merged_Symbol"]
# description
data$Description = id_map[data$Ensembl, "Description"]
# delta pct
data$pct_delta = data$pct.1 - data$pct.2
# ident
data$Ident = "Adipocyte"


### 3. gene quality control
# pct_delta p_val_adj avg_log2FC
data = subset(data, pct_delta>0.1 & p_val_adj<0.05 & avg_log2FC>0.5)


### 4. point target cell type
# all cell types
celltype = levels(ob.merge.rmdb$celltype_final) # 全部细胞类型
# number of cell types
Nc = length(celltype)
# target cell type
target = c("Adipocyte") 
# other cell types
others = celltype[!celltype %in% target] 


### 5. extract average gene expression based on cell types
# cell types - gene expression matrix
expression = AverageExpression(
  ob.merge.rmdb, group.by = 'celltype_final', features = row.names(data),
  assays = "RNA", slot = "data") %>% as.data.frame()
colnames(expression) = celltype
# target expression
exp_target = expression[, target] %>% as.numeric()
# others expression
exp_others = expression[, others]


### 6. calculate specificity index (Si)
# 计算特异性指数 Si
exp_others$Si = NULL
Si = numeric()
for (k in 1:nrow(expression)){
  Si[k] = Nc - sum(exp_others[k,]>0.3*exp_target[k])
  # 注意这里直接将Si 值赋给 express_others$Si[k] 会出现 bug，造成每个 Si 值 +1 
}
exp_others$Si = Si


### 7. grade classification
exp_others$grade = NULL
for (k in 1:nrow(expression)){
  if(exp_others$Si[k] == Nc){
    exp_others$grade[k] = c("A+")
  } else if (exp_others$Si[k] >= (Nc - 0.1*Nc) & exp_others$Si[k] < Nc) {
    exp_others$grade[k] = c("A")
  } else if (exp_others$Si[k] >= (Nc - 0.3*Nc) & exp_others$Si[k] < (Nc - 0.1*Nc)){
    exp_others$grade[k] = c("A-")
  } else if (exp_others$Si[k] < (Nc - 0.3*Nc)){
    exp_others$grade[k] = c("N")
  }
}
# 这种方法的缺点是没有把 lineage 组内差异考虑进来
# 还有就是没有考虑到组内基因表达方差变异和细胞数量


### 8. merge data
head(row.names(exp_others))
head(row.names(data))
data = data.frame(row.names = row.names(data),
                  data, exp_others, target = exp_target)


### 9. final markers
data = subset(data, grade %in% c("A+", "A", "A-"))
data = arrange(data, -Si,-avg_log2FC)
table(data$grade)
marker_adipocyte = data


