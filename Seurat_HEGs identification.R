library(Seurat)
library(zoo)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(AnnotationHub)
library(AnnotationDbi)
library(ComplexHeatmap)
library(circlize) # colorRamp2

### 1. highly expressed genes identification
# *** logfc.threshold = 0.5 ***
all_DEGs_celltype = FindAllMarkers(
  ob.merge.rmdb, logfc.threshold = 0.5, min.pct = 0, only.pos = T)
all_DEGs_celltype$gene = row.names(all_DEGs_celltype)
# transfer to human homologous genes
# self-build sheep-human id transformation file
all_DEGs_celltype$ensembl=id_map_seurat[all_DEGs_celltype$gene,"Ensembl"]
# 因为不同细胞类型筛选到相同基因，例如几个 psc，gene symbo 被重新编码，所以 ensembl 中被填充为 na
# 不能直接删除 na，因为后边要做富集分析
# 先排序，然后用同一列，上一个非 na 值填充
all_DEGs_celltype = arrange(all_DEGs_celltype, gene)
all_DEGs_celltype <- na.locf(all_DEGs_celltype)
all_DEGs_celltype$gene = id_map[all_DEGs_celltype$ensembl,"merged_Symbol"]
all_DEGs_celltype$description = id_map[all_DEGs_celltype$ensembl, "Description"]
# delta pct
all_DEGs_celltype$pct_delta = all_DEGs_celltype$pct.1 - all_DEGs_celltype$pct.2
# p_val_adj
all_DEGs_celltype = subset(all_DEGs_celltype, p_val_adj < 0.05)
all_DEGs_celltype = arrange(all_DEGs_celltype, cluster)


### 2. Functional enrichment analysis
# human database
library(org.Hs.eg.db)
#detach("package:org.Hs.eg.db",unload = TRUE)
# orgdb = loadDb(file = "org.Oar.eg.db") # 读取

# save enrichment results to the all_DEGs_celltype_go dataframe
all_DEGs_celltype_go = data.frame()
obj.markers=all_DEGs_celltype
# 生成n个list，并重命名
gb <- obj.markers %>% group_by(cluster)
l <- group_split(gb)
k <- group_keys(gb)
names(l) <- k[['cluster']]
for (x in names(l)){
  d0 = l[[x]]
  # GO
  tryCatch(
    { #BP
      ego_BP <- enrichGO(gene=d0$gene,
                         OrgDb = org.Hs.eg.db, # or org.Oar.eg.db
                         keyType='SYMBOL', #ENTREZID, ENSEMBL
                         ont = "BP",
                         pAdjustMethod = 'BH',
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05)
      ego_BP <- simplify(ego_BP) %>% as_tibble() %>% mutate(Links=str_c('http://amigo.geneontology.org/amigo/term/',ID))
      ego_BP$ONTOLOGY <- "BP"
      #MF
      ego_MF <- enrichGO(gene=d0$gene,
                         OrgDb = org.Hs.eg.db,
                         keyType='SYMBOL',
                         ont = "MF",
                         pAdjustMethod = 'BH',
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05)
      ego_MF <- simplify(ego_MF) %>% as_tibble() %>% mutate(Links=str_c('http://amigo.geneontology.org/amigo/term/',ID))
      ego_MF$ONTOLOGY <- "MF"
      #CC
      ego_CC <- enrichGO(gene=d0$gene,
                         OrgDb = org.Hs.eg.db,
                         keyType='SYMBOL',
                         ont = "CC",
                         pAdjustMethod = 'BH',
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05)
      ego_CC <- simplify(ego_CC) %>% as_tibble() %>% mutate(Links=str_c('http://amigo.geneontology.org/amigo/term/',ID))
      ego_CC$ONTOLOGY <- "CC"
      #merge
      ego <- rbind(ego_BP, ego_MF, ego_CC)
      ego$cluster = unique(d0$cluster) %>% as.character()
      ego = arrange(ego, qvalue)
      all_DEGs_celltype_go = rbind(all_DEGs_celltype_go, ego)
      #stop run, if no go
      if (dim(ego)[1]==0){
        notice = paste0('No GO in ', x)
        print(notice)
      } else {
        write_tsv(ego, paste0('./degs_celltype_enrichment/',x, '_cluster', '_go.xls'))
        #plot
        ego$LOG = -log10(ego$qvalue)
        ego = ego[order(ego$ONTOLOGY, -ego$LOG),]
        ego_top10 <- ego %>% group_by(ONTOLOGY) %>% top_n(n = 10, wt = LOG)
        #首先需要将go term转换为因子
        ego_top10$Description=factor(ego_top10$Description,levels=ego_top10$Description)
        #重新排序
        ego_top10=ego_top10[order(ego_top10$LOG),]
        p=ggplot(ego_top10, mapping = aes(x=Description, y=LOG, fill=ONTOLOGY))
        #绘制条形
        p0=p+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(ego_top10$Description)))
        ggsave(paste0('./degs_celltype_enrichment/', x, '_cluster', '_go_barplot.png'), 
               width = 7, height = 4, p0, dpi=300)
      }
    }
  )
}
# save
write.csv(all_DEGs_celltype, './degs_celltype/all_DEGs_celltype.csv')
write.csv(all_DEGs_celltype_go, './degs_celltype/all_DEGs_celltype_go.csv')


### 3. generate heatmap using ComplexHeatmap package
# based on all_DEGs_celltype, average gene expression matrix and representively enriched GO terms
library(ComplexHeatmap)
library(circlize) # colorRamp2

# 1. degs
all_DEGs_celltype = read.csv('./degs_celltype/all_DEGs_celltype.csv')
markers = unique(all_DEGs_celltype$ensembl)
id_map_seurat = readRDS('./gene_annotation/id_map_seurat.rds')
row.names(id_map_seurat) = id_map_seurat$Ensembl
markers = id_map_seurat[markers,]

# 2. average gene expression matrix
mean_gene_exp <- as.matrix(
  data.frame(Seurat::AverageExpression(ob.merge.rmdb, 
                                       features = markers$Symbol_uniq, 
                                       group.by = 'celltype_final', 
                                       assays = "RNA", slot = "data")))
colnames(mean_gene_exp) <- levels(ob.merge.rmdb$celltype_final)
# data normalization
data <- t(scale(t(mean_gene_exp), scale = T, center = T))
# cell * gene 矩阵
data = t(data)

# 3. define colors
col_fun = colorRamp2(c(-2, 0, 3), c("#2c1dff", "white", "#ff2d00"))
# 细胞颜色 keys = valuses 格式
# 颜色和数据一定要严格对应
cols_final = c(
  "PSC_C1"="#b1ff91", "PSC_C2"="#2e7d32", "PSC_C3"="#a9cf54", 
  "PSC_C4"="#66bb6a", "PSC_C5"="#43a047", "PSC_C6"="#96ed89",
  "CTP_C1"="#ffc682", "CTP_C2"="#fbc9c9","CTP_C3"="#f57777",         
  "Preadipocyte"="#d23600", "Adipocyte"="#B22222",
  "VSMC_C1"="#ffbe00", "VSMC_C2"="#fff176", 
  "Preosteoblast"="#b8a3de","Osteoblast"="#8a23cd",
  "Prechondrocyte"="#edd4fe", "Chondrocyte_C1"="#dba9fd", "Chondrocyte_C2"="#43026f", 
  "PSC_Myo"="#add5f7", "Satellite_Cell"="#799ae0", "Myoblast"="#1c3ffd", "Myocyte"="#020873",
  "Unknow"="#1c1d21", "Endothelium"="#435862", "Epithelium"="#5d5100",         
  "Erythroblast"="#d9d9d9", "T_Cell"="#dae2e5", "Macrophage"="#99aeb8", 
  "Mast_Cell"="#4a606b", "Neuron"= "#04bfbf")

# 4. left annotation
# cell label and color modulr
celltype = levels(ob.merge.rmdb$celltype_final) 
celltype
ha_left = rowAnnotation(
  # 单独指定细胞label及其对应的颜色模块
  # 1. label
  cell_text = anno_text(celltype, location=1, 
                        just='right', gp = gpar(fontsize = 10)), # 右对齐，字体大小
  # 2. 单独添加颜色快，用 label 映射颜色快
  Cell_type = celltype,
  col = list(Cell_type = cols_final),
  # 是否显示 cell_text 和 Cell_type 两种注释信息图例
  show_legend = c(F, F)) 
ha_left

# 5. right annotation
# representively enriched GO termsead.csv('./degs_celltype/all_DEGs_celltype_go_plot.csv')
go$p_adjust = -log10(go$p.adjust) # 矫正 p_adjust
row.names(go) = go$cluster
go = go[, c('Description', 'p_adjust')]
# go label
go_text = data.frame(row.names = row.names(go), Description = go$Description)
# go p_adjust
go_value = data.frame(row.names = row.names(go), Description = go$p_adjust)
# 注释
ha_right = rowAnnotation(
  # barplot 注释 go_value
  # fill 填充颜色，width 宽度
  p_adjust = anno_barplot(go_value, gp = gpar(fill = cols_final), width = unit(2, "cm")),
  # text 注释 go label
  go_text = anno_text(go_text$Description,just = 'left', gp = gpar(fontsize = 8)),
  show_legend = c(F, F)) # 去掉图例
ha_right

# 6. plot
pdf('./figures/degs_celltype_heatmap_go.pdf', width = 10, height = 11)
Heatmap(data, col = col_fun, width = unit(6, "cm"), height = unit(11, "cm"),
        # heatmap_width = unit(10, "cm"), heatmap_height = unit(1, "cm"), # heatmap 所有组分整体高度宽度
        cluster_rows = F, show_row_names = F, cluster_columns = F, show_column_names = F,
        left_annotation = ha_left, # 左侧注释
        right_annotation = ha_right, # 右侧注释
        raster_quality = 15, # raster 质量
        show_heatmap_legend = T,
        heatmap_legend_param = list( # legend 的一些参数，但无法设定 legend 位置
          title = c('Exp'), title_position='leftcenter',
          # topleft, topcenter, leftcenter-rot and lefttop-rot for vertical legend
          # leftcenter, lefttop are only for horizontal legend
          legend_direction = c('horizontal') # "vertical", "horizontal
        )
        )
dev.off()
