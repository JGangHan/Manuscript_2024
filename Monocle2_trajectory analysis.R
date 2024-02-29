############################# differentiation trajectory analysis ##########################
# cell types: myogenic progenitors (PSC_Myo), satellite cells, myoblasts, myocytes
# color configuration: "PSC_Myo"="#add5f7", "Satellite_Cell"="#799ae0", "Myoblast"="#1c3ffd", "Myocyte"="#020873",
# note: should perform this analysis in server rather than run these codes in local computer, due to long run time and memory

library(Seurat)
library(monocle)
library(dplyr)

### 1. data subset
myo.cell=c("PSC_Myo", "Satellite_Cell", "Myoblast", "Myocyte")
myo.cell= subset(ob.merge.rmdb, celltype_final %in% myo.cell) %>% FindVariableFeatures()
Idents(myo.cell)


### 2. create monocle object from Seurat object
assay="RNA"
pd <- new("AnnotatedDataFrame", data = obj@meta.data)
gene_annotation=data.frame(gene_short_name = rownames(obj[[assay]]),stringsAsFactors=F)
rownames(gene_annotation)<-gene_annotation$gene_short_name
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(GetAssayData(obj,slot="counts",assay=assay),
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())


### 3. construct differentiation trajectory
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
pData(HSMM)$celltype_final = droplevels(pData(HSMM)$celltype_final)
# set ordering genes
HSMM<-setOrderingFilter(HSMM,ordering_genes = VariableFeatures(obj))
# dimension reduction
HSMM <- reduceDimension(HSMM, max_components = 2,reduction_method = 'DDRTree', verbose = T)
# order cells according to ordering genes
HSMM <- orderCells(HSMM)
# opt: re-order cells
HSMM = orderCells(HSMM, root_state = 2)


### 4. plot, overall trajectory
pdf("trajectory_myo_overall.pdf", height =3, width = 2.6)
# state
plot_cell_trajectory(HSMM, color_by = "State",cell_size = 0.5)
# time
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 0.5)
# celltype
cols = c("PSC_Myo"="#add5f7", "Satellite_Cell"="#799ae0", 
         "Myoblast"="#1c3ffd", "Myocyte"="#020873")
plot_cell_trajectory(HSMM, color_by = "celltype_final",cell_size = 0.5) + 
  scale_color_manual(values=cols)
# sample
cols_stage_umap = c("#e5fee7", "#d1eed3", "#bedec0", '#aacfad', '#97bf9a', "#84b088","#72a176", "#5f9264", '#4d8353')
plot_cell_trajectory(HSMM, color_by = "samplename",cell_size = 0.5) + 
  scale_color_manual(values=cols_stage_umap)
dev.off()


### 5. plot, seperated trajectory
# celltype
pdf("trajectory_myo_seperated_celltype.pdf", width = 4, height = 4.5)
plot_cell_trajectory(HSMM, color_by = "celltype_final",cell_size = 1) + 
  scale_color_manual(values=cols) + facet_wrap(~celltype_final, ncol=2)
dev.off()
# sample
pdf("trajectory_myo_seperated_sample.pdf", width = 5.5, height = 6)
plot_cell_trajectory(HSMM, color_by = "samplename",cell_size = 1) + 
  scale_color_manual(values=cols_stage_umap) +
  facet_wrap(~samplename, ncol=3)
dev.off()


### 6. save data
meta.data<-pData(HSMM)[, c("Size_Factor", "Pseudotime", "State")]
obj<-AddMetaData(obj,metadata = meta.data)
saveRDS(obj,"myo_seurat.rds")
saveRDS(HSMM,"myo_monocle2.rds")


### 7. pseudotime-dependent DEGs
diff_pseudotime <- differentialGeneTest(HSMM, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_pseudotime = subset(diff_pseudotime, qval<1.0E-5)
diff_pseudotime = arrange(diff_pseudotime, qval)
write.csv(diff_pseudotime, "./degs_myo_pseudotime.csv")
# gene annotation
data = degs_myo_pseudotime
data$Gene = row.names(data)
# ensembl
data$Ensembl=id_map_seurat[data$Gene,"Ensembl"]
# gene_merged
data$Gene_merged = id_map[data$Ensembl, "merged_Symbol"]
# description
data$Description = id_map[data$Ensembl, "Description"]
# ident
data$Ident = "degs_myo_pseudotome"
degs_myo_pseudotime = data
write.csv(degs_myo_pseudotime, './monocle/degs_myo_pseudotime.csv')


### 8. branch dependent genes (bdGenes)
# state1: myocyte branch
# state3: satellite cell branch
BEAM_res <-BEAM(HSMM, branch_point = 1, branch_states = c(1, 3), cores = 50)
# gene quality control
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(
  subset(fData(HSMM),num_cells_expressed > 0.1*min(table(obj$celltype_final))))
BEAM_res = BEAM_res[row.names(BEAM_res) %in% expressed_genes, ]


### 9. the distance of genes among two branches
# state1: myocyte branch
# state3: satellite cell branch
# 必须要轨迹中已有的branch labels
HSMM_abc = calABCs(HSMM[expressed_genes,], branch_point = 1, 
                   trajectory_states = c(1, 3), # branch_labels = c(State1, State3),
                   cores = 50)


### 10. merge HSMM_abc and BEAM_res
BEAM_res = merge(BEAM_res, HSMM_abc[, c("ABCs", "gene_short_name")], 
                 by = "gene_short_name")
row.names(BEAM_res) = BEAM_res$gene_short_name
BEAM_res = BEAM_res[,c('gene_short_name', 'pval', 'qval', 'ABCs')]
head(BEAM_res)


### 11. quality control
BEAM_res = subset(BEAM_res, qval<1.0E-5 & abs(ABCs) > 5) 
dim(BEAM_res)


### 13. heatmap for all bdGenes
# state 1: Myocyte; state 3: Satellite cell
branched_heatmap = plot_genes_branched_heatmap(
  HSMM[row.names(BEAM_res),], branch_point = 1, branch_states =  c(1, 3),
  branch_labels = c("Myocyte", "Satellite cell"), branch_colors = c("#add5f7", "#020873", "#799ae0"),
  num_clusters = 4,###对基因进行聚类
  cores = 50, show_rownames = F, use_gene_short_name = T, return_heatmap=T)
ggsave(filename = 'branched_heatmap_myo.pdf',
       plot = branched_heatmap$ph_res, 
       width = 4, height = 4, device = cairo_pdf)

# 热图 cluster 信息补充到 BEAM_res
cluster = branched_heatmap$annotation_row
table(cluster$Cluster)
for (i in 1:nrow(BEAM_res)){
  BEAM_res[i,"cluster"] = cluster[row.names(BEAM_res[i,]),1]
}
head(BEAM_res)
write.csv(BEAM_res, "bdgs_myo.csv")


### 14. gene annotation of bdGenes
bdgs_myo = read.csv('./monocle/bdgs_myo.csv', row.names = 1)
data = bdgs_myo
# gene
data$Gene = row.names(data)
# ensembl
data$Ensembl=id_map_seurat[data$Gene,"Ensembl"]
# gene_merged
data$Gene_merged = id_map[data$Ensembl, "merged_Symbol"]
# description
data$Description = id_map[data$Ensembl, "Description"]
# ident
data$Ident = "myo branched dependent genes"
bdgs_myo = data
bdgs_myo = arrange(bdgs_myo, qval)




