setwd('./pre_study')

#BiocManager::install('Seurat',update=F)
library(Seurat)
gc()
# 读取数据
gc()
#读取数据
Normal <- Read10X(data.dir = './Normal/')
Tumor <- Read10X(data.dir = './Tumor/')

#分别创建SeuratObject
scRNA1 <- CreateSeuratObject(counts = Normal, project = "Normal")
scRNA2 <- CreateSeuratObject(counts = Tumor, project = "Tumor")

#利用merge函数把两个样本合在一起，并且给每个样本都加上各自的样本名，方便后续分析
sce.all = merge(scRNA1, y = c(scRNA2), add.cell.ids = c("Normal", "Tumor"),
                project = 'ESCA', merge.data = TRUE)
scRNA=sce.all
# metadata为样本信息，我们需要定义分组
scRNA@meta.data$tissue_type=scRNA@meta.data$orig.ident
#install.packages('stringr')
# 去除数字,获得干净的分组
scRNA@meta.data$tissue_type=stringr::str_remove(scRNA@meta.data$tissue_type,'[0-9]')

# 质控
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
library(ggplot2)
col.num <- length(levels(scRNA@active.ident))
# 过滤前
VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 


##设置质控标准，比较随意
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=4000
pctMT=15

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#降维聚类#########################
library(Seurat)
library(tidyverse)
library(patchwork)

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
plot

#只对高变基因进行scale
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)


scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
plotc


#PC选取比较随意
pc.num=1:15

scRNA <- FindNeighbors(scRNA, dims = pc.num) 

scRNA <- FindClusters(scRNA)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

#tSNE可视化
scRNA = RunTSNE(scRNA, dims = pc.num)
plot1 = DimPlot(scRNA, reduction = "tsne",label = T, group.by="orig.ident") 
plot1

#UMAP可视化
scRNA <- RunUMAP(scRNA, dims = pc.num)
plot2 = DimPlot(scRNA, reduction = "umap",label=T,split.by ="orig.ident" ) 
plot2

#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc



###鉴定细胞类型

habermann_imm <- c('CD274',"CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3",'GATA3', "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")

DotPlot(scRNA, features = habermann_imm,group.by = "seurat_clusters") + coord_flip()
DotPlot(scRNA, features = habermann_oth,group.by = "seurat_clusters") + coord_flip()

### 注释

#0-3
#4-7
#8-11
#12-15
#16-18
celltype <- c('B_cells' ,'NK_cells','T_cells','Plasma_cells',
              'T_cells','Monocytes','Macrophages', 'T_cells ', 
              'Macrophages','B_cells','Epithelial' ,'NK_cells',
              'T_cells','Fibroblasts','B_cells','Macrophages', 
              'Mast_cells','B_cells','NK_cells')


#celltype <- c('B_cells','T_cells','Epithelial_cells','T_cells','Endothelial_cells',
#              'B_cells','Fibroblasts','Endothelial_cells','Myeloids',
#              'Epithelial_cells','B_cells','Epithelial_cells','Epithelial_cells')

Idents(scRNA) <- scRNA@meta.data$seurat_clusters
names(celltype) <- levels(scRNA)
scRNA<- RenameIdents(scRNA, celltype)

scRNA@meta.data$celltype <- Idents(scRNA)
#!!!!!!
Idents(scRNA)=scRNA@meta.data$celltype
#library("ggsci")
#mypal=pal_npg("nrc",alpha=0.7)(10)
#library("scales")
#mypal
colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4', "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2")
p1 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='tsne',pt.size = 1,cols = colors,split.by ="orig.ident")
p1
p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='umap',pt.size = 1,cols = colors,split.by="orig.ident")
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
p3
FeaturePlot(scRNA,features = 'TUBA1A',split.by = 'tissue_type',label=T)
saveRDS(scRNA,file ='scRNA_anno.rds')

#scRNA=scRNA_anno
#scRNA$celltype=gsub('T_cells ','T_cells',scRNA$celltype)
#scRNA$celltype=gsub('Macrophages ','Macrophages',scRNA$celltype)
#table(scRNA$celltype)
## 一定读取我们准备的！！！！
scRNA=readRDS('scRNA_anno.rds')


#########细胞通讯全局分析
# 先安装cellchat,4.1版本 的R为宜
#devtools::::install_github("sqjin/CellChat")
library(CellChat)
# 选取第一个样本，如果内存不够，再进一步减少细胞数，例如随机抽1000个
# 内存够，则不挑选直接上
#scRNA_chat <- subset(scRNA, orig.ident=='Tumor1')
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

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

write.csv(df.net,file ='cellchat.csv',quote=F)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()
??netVisual_circle
#####------------------------
# 关键热图

library(Seurat)
gene=read.table('../REACTOME_AGGREPHAGY.v7.5.1.gmt')
gene=gene[,3:ncol(gene)]
gene=t(gene)
DoHeatmap(subset(scRNA,downsample=50,),features = gene,group.by = 'celltype',assay='RNA',slot = 'data',
          group.colors =c('#313c63','#b42e20','#ebc03e','#377b4c',
                     '#7bc7cd','#5d84a4'),lines.width = 10)+
  scale_fill_gradientn(colors=c('white','firebrick3'),na.value = 'white')
metadata=scRNA@meta.data


VlnPlot(scRNA,features = c('HSP90AA1','	RPS27A','UBA52','UBB','UBC',"VIM"))
VlnPlot(scRNA,features = 'HSP90AA1')
FeaturePlot(scRNA,features = 'HSP90AA1',reduction = 'umap',pt.size = 1)


