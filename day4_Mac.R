setwd('xxxxx')
load('./scRNA_anno.rds')
scRNA=scRNA_anno
library(Seurat)
DimPlot(scRNA,split.by = 'tissue_type',reduction = 'tsne')
scRNA$celltype=gsub('Macrophages ','Macrophages',scRNA$celltype)
# 巨噬细胞的NMF和分析
scRNA_myelo=subset(scRNA,celltype %in% 'Macrophages')

dim(scRNA_myelo)

## 重新聚类分群
scRNA_myelo<- ScaleData(scRNA_myelo)

scRNA_myelo<- RunPCA(scRNA_myelo) 
plot1 <- DimPlot(scRNA_myelo, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA_myelo, ndims=20, reduction="pca") 
plotc <- plot1+plot2
plotc

#PC选取比较随意
pc.num=1:15

scRNA_myelo <- FindNeighbors(scRNA_myelo, dims = pc.num) 
scRNA_myelo <- FindClusters(scRNA_myelo)
table(scRNA_myelo@meta.data$seurat_clusters)

#tSNE可视化
scRNA_myelo = RunTSNE(scRNA_myelo, dims = pc.num)
plot1 <-DimPlot(scRNA_myelo, reduction = "tsne",label = T,split.by = 'orig.ident') 

#UMAP可视化
scRNA_myelo <- RunUMAP(scRNA_myelo, dims = pc.num)
plot2 <-DimPlot(scRNA_myelo, reduction = "umap",label=T,split.by = 'orig.ident') 
plotc <- plot1+plot2
plotc
## 髓系marker
habermann_imm <- c( "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C",
                    "CPA3",'GATA3', "KIT", "MKI67", "CDK1")
## 可以看到1，2是巨噬细胞，0是肥大细胞，3是单核细胞
library(ggplot2)
DotPlot(scRNA_myelo, features = habermann_imm,group.by = "seurat_clusters") +
  coord_flip()+scale_colour_gradientn(colours =c('#dadada','#bc3c29'))

## 
Myeloid_subtype <- c('Mac','Mac','Mac','Mac','Mono','Mono','Mono','Mono')

Idents(scRNA_myelo) <- scRNA_myelo@meta.data$seurat_clusters
names(Myeloid_subtype) <- levels(scRNA_myelo)
scRNA_myelo<- RenameIdents(scRNA_myelo, Myeloid_subtype)

scRNA_myelo@meta.data$Myeloid_subtype <- Idents(scRNA_myelo)

Idents(scRNA_myelo)=scRNA_myelo@meta.data$Myeloid_subtype

DimPlot(scRNA_myelo, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                  '#7bc7cd','#5d84a4','#dadada'),
        reduction = 'tsne',split.by = 'orig.ident') 



## 取出肿瘤
scRNA_myelo=subset(scRNA_myelo,orig.ident %in% 'Tumor')
## 取出mac
scRNA_mac=subset(scRNA_myelo,Myeloid_subtype %in% 'Mac')


gene=read.table('./REACTOME_AGGREPHAGY.v7.5.1.gmt')
gene=gene[,3:ncol(gene)]
gene=t(gene)

# 聚集体自噬aa
scRNA_mac$id=colnames(scRNA_mac)
scRNA_mac_aa=scRNA_mac[rownames(scRNA_mac) %in% gene,]

df <- scRNA_mac_aa@assays$RNA@data
df=as.data.frame(df)

df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]

## 一点没表达的细胞去除
scRNA_mac=subset(scRNA_mac,id %in% colnames(df))

# nmf
library(NMF)
res <- nmf(df, 10, method = "snmf/r", seed = 'nndsvd')

## NMF结果返回suerat对象
scRNA_mac@reductions$nmf <- scRNA_mac@reductions$pca
scRNA_mac@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_mac@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类!重要
set.seed(999)

scRNA_mac.nmf <- RunUMAP(scRNA_mac, reduction = 'nmf', dims = 1:10) 
scRNA_mac.nmf <- FindNeighbors(scRNA_mac.nmf,reduction = 'nmf', dims = 1:10)
scRNA_mac.nmf <- FindClusters(scRNA_mac.nmf)
scRNA_mac.nmf$NMF_cluster=scRNA_mac.nmf$seurat_clusters

## 结果可视化  
DimPlot(scRNA_mac.nmf, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                      '#7bc7cd','#5d84a4','#dadada')) 

df=FindAllMarkers(scRNA_mac.nmf,logfc.threshold = 0.5,only.pos = T)
sss=intersect(df$gene,gene)
write.csv(df,file ='deg_mac.csv',quote=F)

# test some genes
FeaturePlot(scRNA_mac.nmf,features = 'VCP',label = T)
DoHeatmap(scRNA_mac.nmf,features = rownames(df),assay = 'RNA',slot = 'count')


## 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install()
#BiocManager::install('monocle',update = F,ask = F)

library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNA_mac.nmf@assays$RNA@counts)

data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_mac.nmf@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)

ss=intersect(gene,rownames(scRNA_mac.nmf))
dev.off()
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ss,],
                                                 # num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_pseudotime_cluster 


#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1

plot4 <- plot_cell_trajectory(mycds, color_by = "NMF_cluster")
plot4

##合并出图
plotc <- plot1|plot4
plotc


## 定义新的细胞亚群！注意修改代码！！！！！！！！！！！！！！
NMF_celltype <- c('TUBA1B+Mac-C0','Non-Aggre-Mac-C1','Non-Aggre-Mac-C1','UBB+Mac-C2')

Idents(scRNA_mac.nmf) <- scRNA_mac.nmf@meta.data$NMF_cluster
names(NMF_celltype) <- levels(scRNA_mac.nmf)
scRNA_mac.nmf<- RenameIdents(scRNA_mac.nmf, NMF_celltype)

scRNA_mac.nmf@meta.data$NMF_celltype <- Idents(scRNA_mac.nmf)

Idents(scRNA_mac.nmf)=scRNA_mac.nmf@meta.data$NMF_celltype

DimPlot(scRNA_mac.nmf,group.by = 'NMF_celltype',
        cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
               '#7bc7cd','#5d84a4'),
        label = T,label.size = 3)

saveRDS(scRNA_mac.nmf,file ='scRNA_mac_NMF.RDS')

########## 细胞通讯(巨噬细胞和上皮）#################################
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')
scRNA_chat <-subset(scRNA_chat, celltype=='Epithelial')

# 赋一列
scRNA_chat$NMF_celltype='Epithelial_cells'
scRNA_chat=merge(scRNA_chat,scRNA_mac.nmf)

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "NMF_celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T,sources.use = 'Epithelial_cells',
                 title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, 
                 title.name = "Interaction weights/strength")


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2


#### scMetabolism评估巨噬细胞代谢活性
#devtools::install_github("YosefLab/VISION")
#devtools::install_github("wu-yc/scMetabolism")

library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_mac.nmf_meta<-sc.metabolism.Seurat(obj = scRNA_mac.nmf, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

input.pathway <- rownames(scRNA_mac.nmf_meta@assays[["METABOLISM"]][["score"]])[1:30]
DotPlot.metabolism(obj = scRNA_mac.nmf_meta,
                   pathway = input.pathway, phenotype = "NMF_celltype", norm = "y")
                   

###M1,M2亚型比较！
scRNA_mac.nmf=readRDS('./scRNA_mac_NMF.RDS')
library(dplyr)
library(readr)
library(stringr)
m1m2_pws <- read_lines("CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

### M2
m1m2_pws <- append(m1m2_pws, read_lines("CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.gmt") %>%
                     lapply(str_split, "\\t") %>% 
                     unlist(recursive = F) %>% 
                     lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
                     unlist(recursive = F))

save(m1m2_pws,file ='m1m2.Rdata')
load('m1m2.Rdata')

library(Seurat)
scRNA_mac.nmf<- AddModuleScore(object = scRNA_mac.nmf, 
                               features = m1m2_pws, name = c("m1up", "m1dn"), nbin = 12)


VlnPlot(scRNA_mac.nmf, features = c("m1up1", "m1dn2"), 
        group.by = "NMF_celltype", cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                          '#7bc7cd','#5d84a4'))

library(viridis)
FeaturePlot(scRNA_mac.nmf,features = c("m1up1",'m1dn2'),label = T,cols = magma(10))
