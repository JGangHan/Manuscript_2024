### 1. overall-circular network
groupSize <- as.numeric(table(cellchat@idents))
color.use = c(
  "PSC_C1"="#b1ff91", "PSC_C2"="#2e7d32", "PSC_C3"="#a9cf54", 
  "PSC_C4"="#66bb6a", "PSC_C5"="#43a047", "PSC_C6"="#96ed89",
  "CTP_C1"="#ffc682", "CTP_C2"="#fbc9c9","CTP_C3"="#f57777",         
  "Preadipocyte"="#d23600", "Adipocyte"="#B22222",
  "VSMC_C1"="#ffbe00", "VSMC_C2"="#fff176", 
  "Preosteoblast"="#b8a3de","Osteoblast"="#8a23cd",
  "Prechondrocyte"="#edd4fe", "Chondrocyte_C1"="#dba9fd", "Chondrocyte_C2"="#43026f", 
  "PSC_Myo"="#add5f7", "Satellite_Cell"="#799ae0", "Myoblast"="#1c3ffd", "Myocyte"="#020873",
  "Endothelium"="#435862", "Epithelium"="#5d5100",         
  "Erythroblast"="#d9d9d9", "T_Cell"="#dae2e5", "Macrophage"="#99aeb8", 
  "Mast_Cell"="#4a606b", "Neuron"= "#04bfbf")

pdf("./figures/cellchat_circular_aggregated interactions.pdf", width = 5, height = 5)
# vertex 指的是代表每个 cell group 的点
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 color.use = color.use, vertex.label.cex = 0.75, 
                 edge.width.max = 5, top = 0.05, label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 color.use = color.use, vertex.label.cex = 0.75, 
                 edge.width.max = 5, top = 0.05, label.edge= F, 
                 title.name = "Interaction weights/strength")
dev.off()


### 2. overall-incoming and outcoming signal heatmap
# 计算 network centrality scores
?netAnalysis_computeCentrality
?netAnalysis_signalingRole_heatmap
cellchat <- 
  netAnalysis_computeCentrality(cellchat, slot.name = "netP", thresh = 0.05)
pdf("./figures/heatmap_outcoming and incoming signals.pdf", width = 12, height = 5)
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat, pattern = "outgoing", color.use = color.use)
ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat, pattern = "incoming", color.use = color.use)
ht1 + ht2
dev.off()


### 3. each enriched pathway-circular network
# 无法描述完全特异的概念，因为存在两种形式，一种是配体细胞，一种是受体细胞
pathways.show.all <- cellchat@netP$pathways
groupSize <- as.numeric(table(cellchat@idents))
for (i in 1:length(pathways.show.all)) {
  # 自动保存作图
  netVisual(cellchat, signaling = pathways.show.all[i], color.use = color.use, 
            vertex.weight = 1, weight.scale = T,
            vertex.label.cex = 0.75, #edge.width.max = 5,
            layout = "circle", out.format = c("png"))
}


### 4. committed adipocyts specifically enriched pathways
levels(cellchat@idents)
# 1. adipo to adipo
p_adipo_to_adipo <- subsetCommunication(
  cellchat, sources.use = c(10:11), targets.use = c(10:11))
table(p_adipo_to_adipo$pathway_name) %>% names()
# "ADIPONECTIN", "LAMININ", "MK", "THBS"   

# 2. others to adipo
p_ot_to_adipo <- subsetCommunication(
  cellchat, sources.use = c(1:9, 12:29), targets.use = c(10:11))
table(p_ot_to_adipo$pathway_name) %>% names()
# "CD46", "LAMININ", "MK", "THBS", "WNT"   

# 3. adipo to others
p_adipo_to_ot <- subsetCommunication(
  cellchat, sources.use = c(10:11), targets.use = c(1:9, 12:29))
table(p_adipo_to_ot$pathway_name) %>% names()
# "ADIPONECTIN", "ANGPTL", "COLLAGEN",  "LAMININ", "MK",  "PTN", "SPP1", "THBS", "VEGF" 


### 5. adipo-related pathways
signaling = c(
  # 1. adipo to adipo
  "ADIPONECTIN", "THBS",   
  # 2. others to adipo
  "CD46", "THBS", "WNT",  
  # 3. adipo to others
  "ADIPONECTIN", "ANGPTL", "PTN", "THBS", "VEGF") %>% unique()
pdf("./figures/cellchat_circular_adipo signaling.pdf", width = 5, height = 5)
# Target Cells: Adipocytes and Preadipocytes
netVisual_aggregate(
  cellchat, signaling = signaling, sources.use = c(10, 11), targets.use = c(1:29),
  edge.width.max = 5,
  top = 1,  color.use = color.use, vertex.label.cex = 0.75, label.edge= F)
# Souce Cells: Adipocytes and Preadipocytes
netVisual_aggregate(
  cellchat, signaling = signaling, sources.use = c(1:29), targets.use = c(10,11),
  edge.width.max = 5,
  top = 1,  color.use = color.use, vertex.label.cex = 0.75, label.edge= F)
dev.off()


### 6. hierarchy plot
# six pathways
signaling = c(
  # 1. adipo to adipo，分析时可以忽略这一部分
  "ADIPONECTIN",
  # 2. others to adipo
  "CD46", "THBS", "WNT",
  # 3. adipo to others
  "ADIPONECTIN", "ANGPTL", "VEGF") %>% unique()
features = c("VEGFA", "FLT1", "KDR",       # VEGF通路，Secreted Signaling
             "THBS3", "THBS4", "CD36",     # THBS通路，ECM-Receptor
             "ANGPTL1", "ITGA1", "ITGB1",  # ANGPTL通路，Secreted Signaling
             'ADIPOQ', 'AdpR2',            # ADIPONECTIN 通路，Secreted Signaling
             'WNT2', 'FZD4', 'LRP6',       # WNT 通路，Secreted Signaling
             'ENSOARG00020021482', 'JAG1'  # CD46 通路，Cell-Cell Contact
)
levels(cellchat@idents)
vertex.receiver = c(1:9, 19)
pdf('./figures/hierarchy_adipo_pathways.pdf', width = 7, height = 7)
for (i in 1:length(signaling)) {
  netVisual_aggregate(
    cellchat, signaling = signaling[i], layout = 'hierarchy', 
    vertex.receiver = vertex.receiver, title.space = 10, edge.width.max = 5,
    top = 1,  color.use = color.use, vertex.label.cex = 0.5, label.edge= F)
}
dev.off()
# LP pairs contribution
for (i in 1:length(signaling)) {
  gg <- netAnalysis_contribution(cellchat, signaling = signaling[i])
  ggsave(filename=paste0(signaling[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


### 7. 6个通路贡献度热图
signaling
# network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


### 8. interactions between committed adipocytes, Endothelium, VSMC and neuron
levels(cellchat@idents)
cellchat@netP$pathways
pdf("./figures/dotplot_cellchat_adipo_to_others.pdf", width = 4.5, height = 2.7)
netVisual_bubble(cellchat, sources.use = c(10, 11), 
                 targets.use = c(4, 8, 9, 10, 11, 12, 13, 23, 29), 
                 signaling = signaling, remove.isolate = T, font.size = 7.5,
                 color.grid = "grey70")
dev.off()
pdf("./figures/dotplot_cellchat_others_to_adipo.pdf", width = 5, height = 2.3)
netVisual_bubble(cellchat, sources.use = c(4, 8, 9, 12, 13, 23, 29), 
                 targets.use = c(10, 11), 
                 signaling = signaling, remove.isolate = T, font.size = 7.5,
                 color.grid = "grey70")
dev.off()


### 9. expression levels of genes encoding LR pairs 
ob.merge.rmdb = readRDS("ob.merge.rmdb.rds")
library(Seurat)
library(ggplot2)
a <- subsetCommunication(
  cellchat, slot.name = "net", signaling = signaling,
  sources.use = c(10, 11), targets.use =  c(4, 8, 9, 10, 11, 12, 13, 23, 29))
b <- subsetCommunication(
  cellchat, slot.name = "net", signaling = signaling,
  sources.use = c(4, 8, 9, 12, 13, 23, 29), targets.use = c(10, 11))
# 所有基因
features = c("VEGFA", "FLT1", "KDR",       # VEGF通路，Secreted Signaling
             "THBS1", "THBS3", "THBS4", "CD36",     # THBS通路，ECM-Receptor
             "ANGPTL1", "ITGA1", "ITGB1",  # ANGPTL通路，Secreted Signaling
             'ADIPOQ', 'AdpR2',            # ADIPONECTIN 通路，Secreted Signaling
             'WNT2', 'FZD4', 'LRP6',       # WNT 通路，Secreted Signaling
             'ENSOARG00020021482', 'JAG1'  # CD46 通路，Cell-Cell Contact
             )

idents = c("Preadipocyte", "Adipocyte", "PSC_C4", "CTP_C2", "CTP_C3", 
           "VSMC_C1", "VSMC_C2", "Endothelium", "Neuron")
cols = c("#66bb6a","#fbc9c9","#f57777","#d23600","#B22222",
         "#ffbe00","#fff176","#435862","#04bfbf")
pdf("./figures/violin_cellchat genes.pdf", width = 5.8, height = 5)
VlnPlot(ob.merge.rmdb, idents = idents, features = features, adjust = 1,
             fill.by = 'ident', cols = cols, 
             stack = T, flip = T) + NoLegend()+
  theme(
    axis.line = element_line(linewidth = 0.5)  # 坐标轴线条粗细
    )
dev.off()
