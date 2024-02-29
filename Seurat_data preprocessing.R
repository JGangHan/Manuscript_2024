library(Seurat)
library(dplyr)


### 1. import dataset from 10X outputs and cell quality control
samplename=c("E50","E55","E60","E63","E66","E69", "E72","E75","E80")
ob.list <- lapply(samplename, function(x) {
  obj = Read10X(paste0(x,"/filtered_feature_bc_matrix"))%>% CreateSeuratObject()
  obj[['percent.mito']] = PercentageFeatureSet(object = obj, pattern = '^(MT|mt|Mt)-')  # mitochondrial gene proportion
  obj$samplename = x
  mt.p <- pnorm(obj$percent.mito, mean = median(obj$percent.mito),  ## 计算中位数
                sd = mad(obj$percent.mito), lower.tail = FALSE)     ## Compute the median absolute deviation
  mt.lim <- min(obj$percent.mito[which(p.adjust(mt.p, method = "fdr") < 0.05)]) # mt threshold
  mito_threshold <-  mt.lim 
  obj = subset(obj, subset = nCount_RNA >= 300 & percent.mito <= mito_threshold) # RNA count number threshold
  obj 
})
percent.mito=ob.merge@meta.data[,c("samplename","percent.mito")]


### 2. gene symbol
names(ob.list)=samplename
id_map=read_delim(paste0(samplename[1] ,'/filtered_feature_bc_matrix/features.tsv.gz'), delim="\t",
                  col_names = c('Ensembl', 'Symbol', 'Type')) %>% 
  dplyr::select(-Type) %>%      ##delete 'Type' colume
  mutate(Symbol_uniq=make.unique(Symbol))


### 3. merge sample without batch effect correction 
ob.merge=merge(ob.list[[1]],ob.list[-1])
table(ob.merge$samplename)


### 4. data pre-processing
ob.merge=NormalizeData(ob.merge)%>%  ## default=LogNormalize
  FindVariableFeatures()%>%          ##genes that highly expressed in some cells, and lowly expressed in others
  ScaleData()%>%                     ####默认的是2000 variable features,全部的features=rownames(ob.merge)
  RunPCA(npcs=100)%>%
  RunUMAP(,dims=1:50)%>%             ##dims must be NULL to run on features
  FindNeighbors(,dims=1:50)%>%       ##features, used only when dims is NULL
  FindClusters()                     ##resolution=0.8


### 5. cell cycle score
sheep_s_gene=c("ATAD2", "BLM", "BRIP1", "CASP8AP2", "CCNE2", "CDC45", "CDC6", "CDCA7")
sheep_g2m_gene=c("ANLN", "ENSOARG00020007614", "AURKA", "AURKB", "ENSOARG00020018892", "BUB1", "CCNB2")
ob.merge=CellCycleScoring(ob.merge,s.features=sheep_s_gene,g2m.features=sheep_g2m_gene)
colnames(ob.merge@meta.data) 


### 6. doublet prediction 
double_rate=c("E50"=0.016,"E55"=0.039,"E60"=0.023,
              "E63"=0.031,"E66"=0.046,"E69"=0.054,
              "E72"=0.054,"E75"=0.054,"E80"=0.061)
samplename=c("E50","E55","E60","E63","E66",
             "E69","E72","E75","E80")
# strategy
# 循环处理9个样本，数据预处理和预测双细胞
# 1) DoubletFinder只能对当个样本预测双细胞，而不能对多个样本整合的单细胞数据预测双细胞
# 2) 去除低质量细胞、低UMI细胞、高下班立体含量细胞
# 3) 进行数据预处理：Normalize、Finvariablefeatures、scale、runpca、runumap等分析
# 4) **对经过预处理的多样本整合数据ob.merge直接拆分成单个样本数据后需要再对单个样本进行预处理
# 5) **主要是因为在多样本中，scale会将基因表达量在所有样本的所有细胞中均值设定为0，
#    其scale后的值当然会与基因表达量在单样本细胞中scale后的值存在出入
names(ob.list)
for (i in names(ob.list)) {
  ob.list[[i]]=NormalizeData(ob.list[[i]])
  ob.list[[i]]=FindVariableFeatures(ob.list[[i]])
  ob.list[[i]]=ScaleData(ob.list[[i]])
  ob.list[[i]]=RunPCA(ob.list[[i]])
  ob.list[[i]]=RunUMAP(ob.list[[i]], dims = 1:10)
}
for(i in names(ob.list) ) {
  sweep.res.list <- paramSweep_v3(ob.list[[i]], PCs = 1:10, sct = F,) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic( ob.list[[i]]$anno)
  DoubletRate=double_rate[i]
  nExp_poi <- round(DoubletRate*ncol(ob.list[[i]])) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ob.list[[i]]  <- doubletFinder_v3(ob.list[[i]] , PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  reuse.pANN = colnames(ob.list[[i]]@meta.data)[dim(ob.list[[i]]@meta.data)[2]-1]
  ob.list[[i]] <- doubletFinder_v3(ob.list[[i]], PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, sct = F)
  
}
doublet=c()
lapply(ob.list,function(x) {
  num_col=ncol(x@meta.data)
  tmp.meta=x@meta.data
  high.doublet = rownames(tmp.meta[tmp.meta[,num_col-1]=="Doublet"&tmp.meta[,num_col]=="Doublet",])
  high.doublet
})%>%do.call("c",.)->doublet  ##合并所有样本的doublet
length(doublet)

# ob.merge中添加doublet属性
ob.merge$doublet="Singlet"
ob.merge@meta.data[doublet, 'doublet']="Doublet"
pdf("umap_doublet_ob.merge.pdf")
DimPlot(ob.merge, group.by = "doublet")
FeaturePlot(ob.merge, features = "nCount_RNA", order=T)  ##默认为灰白色和蓝色"lightgrey and blue"
FeaturePlot(ob.merge, features = "nFeature_RNA", order=T)  ##默认为灰白色和蓝色"lightgrey and blue"
dev.off()
saveRDS(ob.merge,"ob.merge.rds")


### 7. doublet removal and re-data preprocessing
ob.merge.rmdb=subset(ob.merge,subset=doublet!="Doublet")
dim(ob.merge.rmdb)
ob.merge.rmdb = FindVariableFeatures(ob.merge.rmdb)%>%
  ScaleData()%>%
  RunPCA(,npcs=100)%>%
  RunUMAP(,dims=1:50)%>%
  FindNeighbors(,dims=1:50)%>%
  FindClusters()


### 8. plotting
pdf("umap_samplename_clusters_ob.merge.rmdb.pdf",height = 6)
#samplename
DimPlot(ob.merge.rmdb, group.by = "samplename",
        pt.size = 0.1, label = F, label.size = 4,
        label.color = "black", label.box = F, 
        shape.by = NULL, order = NULL, repel = F)##不将细胞群的label与细胞群重叠
#seurat_clusters
DimPlot(ob.merge.rmdb, cols=mycolor, group.by = "seurat_clusters",
        pt.size = 0.1, label = T, label.size = 4,
        label.color = "black", label.box = F, 
        shape.by = NULL, order = NULL, repel = F)##不将细胞群的label与细胞群重叠
dev.off()






