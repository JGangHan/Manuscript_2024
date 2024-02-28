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


