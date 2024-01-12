gc()
load('./scRNA_anno.rds')
scRNA=scRNA_anno
library(Seurat)
DimPlot(scRNA,split.by = 'tissue_type',reduction = 'tsne')
#scRNA$celltype=gsub('T_cells ','T_cells',scRNA$celltype)
# T细胞的NMF和分析
scRNA_T=subset(scRNA,celltype %in% 'T_cells')

dim(scRNA_T)

## T细胞重新聚类分群
scRNA_T <- ScaleData(scRNA_T)

scRNA_T <- RunPCA(scRNA_T) 
plot1 <- DimPlot(scRNA_T, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA_T, ndims=20, reduction="pca") 
plotc <- plot1+plot2
plotc

#PC选取比较随意
pc.num=1:15

scRNA_T <- FindNeighbors(scRNA_T, dims = pc.num) 

scRNA_T <- FindClusters(scRNA_T)
table(scRNA_T@meta.data$seurat_clusters)

#tSNE可视化
scRNA_T = RunTSNE(scRNA_T, dims = pc.num)
plot1 <-DimPlot(scRNA_T, reduction = "tsne",label = T,split.by = 'orig.ident') 

#UMAP可视化
scRNA_T <- RunUMAP(scRNA_T, dims = pc.num)
plot2 <-DimPlot(scRNA_T, reduction = "umap",label=T,split.by = 'orig.ident') 
plotc <- plot1+plot2
plotc
habermann_imm <- c("CD3E", "CD4",  "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7",'IL17A','RORC',"FOXP3")
## 可以看到1，2是巨噬细胞，0是肥大细胞，3是单核细胞
library(ggplot2)
DotPlot(scRNA_T, features = habermann_imm,group.by = "seurat_clusters") + coord_flip()
FeaturePlot(scRNA_T,features = 'CD4',split.by = 'orig.ident')

## 
T_cell_subtype <- c('CD4','CD4','CD4','CD4','CD4',
              'CD8','CD4','CD4','CD8','CD4')

Idents(scRNA_T) <- scRNA_T@meta.data$seurat_clusters
names(T_cell_subtype) <- levels(scRNA_T)
scRNA_T<- RenameIdents(scRNA_T, T_cell_subtype)

scRNA_T@meta.data$T_cell_subtype <- Idents(scRNA_T)

Idents(scRNA_T)=scRNA_T@meta.data$T_cell_subtype

DimPlot(scRNA_T, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                          '#7bc7cd','#5d84a4','#dadada'),
        reduction = 'tsne',split.by = 'orig.ident') 



## 取出肿瘤
scRNA_tumor=subset(scRNA_T,orig.ident %in% 'Tumor')
## 取出CD8
scRNA_T=subset(scRNA_tumor,T_cell_subtype %in% 'CD8')

gene=read.table('./REACTOME_AGGREPHAGY.v7.5.1.gmt')

gene=gene[,3:ncol(gene)]
gene=t(gene)

# 聚集体自噬aa
scRNA_T$id=colnames(scRNA_T)
scRNA_T_aa=scRNA_T[rownames(scRNA_T) %in% gene,]

df <- scRNA_T_aa@assays$RNA@data
df=as.data.frame(df)

df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]

## 一点没表达的细胞去除
scRNA_T=subset(scRNA_T,id %in% colnames(df))

library(NMF)
res <- nmf(df, 10, method = "snmf/r", seed = 'nndsvd')

## NMF结果返回suerat对象
scRNA_T@reductions$nmf <- scRNA_T@reductions$pca
scRNA_T@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_T@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类!重要
set.seed(999)

scRNA_T.nmf <- RunUMAP(scRNA_T, reduction = 'nmf', dims = 1:10) 
scRNA_T.nmf <- FindNeighbors(scRNA_T.nmf,reduction = 'nmf', dims = 1:10)
scRNA_T.nmf <- FindClusters(scRNA_T.nmf)
scRNA_T.nmf$NMF_cluster=scRNA_T.nmf$seurat_clusters

## 结果可视化  
DimPlot(scRNA_T.nmf, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                          '#7bc7cd','#5d84a4','#dadada')) 

df=FindAllMarkers(scRNA_T.nmf,logfc.threshold = 0.5,only.pos = T)



write.csv(df,file ='deg_cd8T.csv',quote=F)

# test some genes
FeaturePlot(scRNA_T.nmf,features = 'TUBB4B',label = T)
DoHeatmap(scRNA_T.nmf,features = rownames(df),assay = 'RNA',slot = 'count')


## 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install()
#BiocManager::install('monocle',update = F,ask = F)

library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNA_T.nmf@assays$RNA@counts)

data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_T.nmf@meta.data)
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

ss=intersect(gene,rownames(scRNA_T.nmf))
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

###
# 若一群细胞TUBB4B top marker,但是倍数0.8，那就是unclear
# 注意合并！！！！！
## 定义新的细胞亚群！注意修改代码！！！！！！！！！！！！！！
NMF_celltype <- c(	'VIM+CD8+T_cells-C0','TUBB4B+CD8+T_cells-C1',
                  'HSP90AA1+CD8+T_cells-C2')

Idents(scRNA_T.nmf) <- scRNA_T.nmf@meta.data$NMF_cluster
names(NMF_celltype) <- levels(scRNA_T.nmf)
scRNA_T.nmf<- RenameIdents(scRNA_T.nmf, NMF_celltype)

scRNA_T.nmf@meta.data$NMF_celltype <- Idents(scRNA_T.nmf)

Idents(scRNA_T.nmf)=scRNA_T.nmf@meta.data$NMF_celltype

DimPlot(scRNA_T.nmf,group.by = 'NMF_celltype',
        cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                                         '#7bc7cd','#5d84a4'),
        label = T,label.size = 3)

saveRDS(scRNA_T.nmf,file ='scRNA_T_NMF.RDS')

########## 细胞通讯(CD8和上皮）#################################
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')
scRNA_chat <-subset(scRNA_chat, celltype=='Epithelial')

# 赋一列
scRNA_chat$NMF_celltype='Epithelial_cells'
scRNA_chat=merge(scRNA_chat,scRNA_T.nmf)

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
                 weight.scale = T, label.edge= T, targets.use = 'Epithelial_cells',
                 title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, 
                 title.name = "Interaction weights/strength")


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2



### 转录因子，记得把前一种细胞的都搬运走,除了保留两个基因组文件外不要放
# 否则覆盖丢失！！！！！
#####SCENIC#####
scRNAsub=readRDS('./scRNA_T_NMF.RDS')
# aertslab/SCENIC
library(SCENIC)
library(Seurat)

## 内存不够控制细胞数目
#set.seed(123456)
#a=sample(1:ncol(scRNAsub),200)
#scRNAsub@reductions$nmf=NULL
#scRNAsub=scRNAsub[,a]

exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo <-scRNAsub@meta.data

# 进入子文件夹，先删去上节课的！！！注意图片保存好！！！！！！！！
setwd("SCENIC") 

Idents(scRNAsub)=scRNAsub$NMF_celltype

### Initialize settings
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

Idents(scRNAsub) <- "NMF_cluster"
# Color to assign to the variables (same format as for NMF::aheatmap)

# 人！！！！
org='hgnc'


dbDir="./SCENIC" # RcisTarget databases location
myDatasetTitle="SCENIC Challenge" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
scenicOptions <- initializeScenic(org=org,
                                  dbDir=dbDir, 
                                  dbs = dbs)

scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"

# 节省内存，我们只取一个基因组文件
# 人：用下面！！！！！！！！！
scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr" = "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v8"
saveRDS(scenicOptions, file="scenicOptions.Rds") 

#Gene filter: Filter by the number of cells in which the gene is detected (minCountsPerGene, by default 6 UMI counts across all samples)
# and by the number of cells in which the gene is detected (minSamples, by default 1% of the cells)
# 尽可能减少基因数目
# 如果需要更多基因，0.1调成0.05或0.01
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.1*ncol(exprMat),minSamples=ncol(exprMat)*.1)
exprMat_filtered <- exprMat[genesKept,]
dim(exprMat_filtered)

#calculating correlation
runCorrelation(exprMat_filtered, scenicOptions)

#GENIE3: infer potential transcription factor targets based on the expression data
exprMat_filtered <- log2(exprMat_filtered+1) 


# 等待较久的一步，纳入基因和细胞越多，等待越久，此处只用了1000多基因，200个细胞
# 对应的TF也只有100多个了，其实有1000多个TF，如果是性能好的计算机可以考虑前面0.1调成0.05或0.01
runGenie3(exprMat_filtered, scenicOptions)


scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123


exprMat_log <- log2(exprMat+1)
dim(exprMat)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# step2容易崩，节省内存只看top50
gc()
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50")) #** Only for toy run!!

# 回到单线程！否则容易报错！！！！！
scenicOptions@settings$nCores <- 1
library(foreach)

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="scenicOptions.Rds") # To save status

#============================================

#下次可以直接读取他,我们已经越过最耗时间的步骤
# 如果经常报错！！！不要重新readRDS和load！！！！应不间断运行

#scenicOptions <- readRDS("scenicOptions.Rds")
#load('scRNAsub.Rdata')
scenicOptions@settings$seed <- 123 # same seed for all of them

#Binarizing the network
runSCENIC_4_aucell_binarize(scenicOptions)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)
colnames(aucell_regulonAUC.t) <- gsub(" ", "_", colnames(aucell_regulonAUC.t))
rownames(aucell_regulonAUC.t) <- gsub("[.]", "-", rownames(aucell_regulonAUC.t))

cd8.scenic <- scRNAsub
# 结果导入endo.scenic
cd8.scenic@meta.data <- cbind(cd8.scenic@meta.data, aucell_regulonAUC.t[rownames(cd8.scenic@meta.data),])
dev.off()
DimPlot(cd8.scenic, reduction = "umap")
# 找6个

#Regulon scores heatmap 热图
Idents(cd8.scenic)=cd8.scenic$NMF_celltype
cells.ord.cluster <- cd8.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)

library(gplots)
library(pheatmap)



#Average Regulon Activity 平均调控活性
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$NMF_celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
library(ggplot2)
library(cowplot)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100))

dev.off()


##### 关键免疫检查点基因的表达
# 回到外面的文件夹!!!!!!!!!!!!!!!!!
setwd('..')
gc()

scRNA=readRDS('./scRNA_T_NMF.RDS')

# 求平均
cd8_avg=AverageExpression(scRNA,group.by = 'NMF_celltype')
cd8_avg=as.data.frame(cd8_avg)

## 根据自己的cluster修改，有些同学可能非常多C
colnames(cd8_avg)=c('C0','C1','C2')
cd8_avg=cd8_avg[,c(1,2,3)]
tme=read.table('co-inhibi.txt')
tme=tme$V1


## ss必须把tme放前面！
ss=intersect(tme,rownames(cd8_avg))
rt=cd8_avg[ss,]
rt=rt[rowMeans(rt) != 0,]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols =F,cluster_rows = F,
                   color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),border_color = NA,
                  )



## CD8 T亚型
library(data.table)
library(readxl)
gene_set = read.table(file='CD8_subtype.txt', header = T, sep = '\t',stringsAsFactors = F)

list <- list()
for(i in 1:length(unique(gene_set$CD8_subtype))){
  list[[i]] <- gene_set$marker[gene_set$CD8_subtype== (unique(gene_set$CD8_subtype)[i])]
}
names(list)<- unique(gene_set$CD8_subtype)

scRNA=readRDS('./scRNA_T_NMF.RDS')
library(Seurat)
names(list)
scRNA=AddModuleScore(scRNA,features = list)

colnames(scRNA@meta.data)
## 注意选择列
df=scRNA@meta.data[,c(13:15)]
a=aggregate(df[,2:3],list(df$NMF_celltype),mean)
rownames(a)=a$Group.1
a$Group.1=NULL
colnames(a)=names(list)
a=t(a)

####自己调整自己的顺序！！！！！！！！！！！！！不一定是5群！！！！
a=a[,c(1,2,3,4,5)]
pheatmap::pheatmap(a,scale = 'row',cluster_cols = F,cluster_rows = F,
                   color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
)
