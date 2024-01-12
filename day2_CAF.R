setwd('e:/work/aggre//')
load('./scRNA_anno.RDS')
library(Seurat)
scRNA=scRNA_anno
# 间质细胞的NMF和分析
# 取出tumor的间质
scRNA_tumor=subset(scRNA,orig.ident %in% 'Tumor')
scRNA_fibro=subset(scRNA_tumor,celltype %in% 'Fibroblasts')
table(scRNA$celltype)
# 查看细胞数目
dim(scRNA_fibro)

gene=read.table('./REACTOME_AGGREPHAGY.v7.5.1.gmt')
gene=gene[,3:ncol(gene)]
gene=t(gene)

# 聚集体自噬aa(用聚集体自噬的表达矩阵！)
scRNA_fibro$id=colnames(scRNA_fibro)
scRNA_fibro_aa=scRNA_fibro[rownames(scRNA_fibro) %in% gene,]

## 抽取表达矩阵
df <- scRNA_fibro_aa@assays$RNA@data
df=as.data.frame(df)

df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]

## 一点没表达的细胞去除
scRNA_fibro=subset(scRNA_fibro,id %in% colnames(df))


# 开始nmf
#install.packages("NMF")  
library(NMF)
res <- nmf(df, 10, method = "snmf/r", seed = 'nndsvd')

## NMF结果返回suerat对象
scRNA_fibro@reductions$nmf <- scRNA_fibro@reductions$pca
scRNA_fibro@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_fibro@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类!重要
set.seed(999)
library(Seurat)
scRNA_fibro.nmf <- RunUMAP(scRNA_fibro, reduction = 'nmf', dims = 1:10) 
scRNA_fibro.nmf <- FindNeighbors(scRNA_fibro.nmf,reduction = 'nmf', dims = 1:10)
scRNA_fibro.nmf <- FindClusters(scRNA_fibro.nmf)

scRNA_fibro.nmf$NMF_cluster=scRNA_fibro.nmf$seurat_clusters

## 结果可视化，群太多自行加颜色
DimPlot(scRNA_fibro.nmf, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                          '#7bc7cd','#5d84a4','#dadada')) 

## 每个cluster的特征基因
df=FindAllMarkers(scRNA_fibro.nmf,logfc.threshold = 0.5,only.pos = T)

write.csv(df,file ='deg_fibro.csv',quote=F)
#NMF_celltype <- c('DYNC1LI2+Tsc-C0','UBB+Tsc-C1')
## 画图看一看
FeaturePlot(scRNA_fibro.nmf,features = 'DYNC1LI2')
FeaturePlot(scRNA_fibro.nmf,features = 'UBB')
# logFC比较大，2.特征基因排在首位的是AA相关基因
#如果无logFC>1,且全AA无关：Non-aggre
# 如果无LOGFC>1,但是AA有关，AA0.7：Unclear

## 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install()
#BiocManager::install('monocle',update = F,ask = F)
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       #'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       #'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       #'terra', 'ggrastr'))
#BiocManager::install("monocle") #作者：小云爱生信 https://www.bilibili.com/read/cv24584412 出处：bilibili
#没有的包先安装
library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNA_fibro.nmf@assays$RNA@counts)

data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_fibro.nmf@meta.data)
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
#排序，报错请用4.1以上的R，并重装monocle
mycds <- orderCells(mycds)

ss=intersect(gene,rownames(scRNA_fibro.nmf))
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


## 定义新的细胞亚群！

## 定义新的细胞亚群！注意修改代码！！！！！！！！！！！！！！
## 把Non-aggre的放到C3
NMF_celltype <- c('TUBA1A+CAF-C0','HSP90AA1+CAF-C1',"VIM+CAF-C2")

#!
Idents(scRNA_fibro.nmf) <- scRNA_fibro.nmf@meta.data$NMF_cluster
names(NMF_celltype) <- levels(scRNA_fibro.nmf)
scRNA_fibro.nmf<- RenameIdents(scRNA_fibro.nmf, NMF_celltype)

scRNA_fibro.nmf@meta.data$NMF_celltype <- Idents(scRNA_fibro.nmf)


Idents(scRNA_fibro.nmf)=scRNA_fibro.nmf@meta.data$NMF_celltype


DimPlot(scRNA_fibro.nmf,group.by = 'NMF_celltype',label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                                           '#7bc7cd','#5d84a4'))

FeaturePlot(scRNA_fibro.nmf,features = 'TUBA1A')
saveRDS(scRNA_fibro.nmf,file ='scRNA_fibro_NMF.RDS')


########## 细胞通讯(间质和上皮）#################################
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')
scRNA_chat <-subset(scRNA_chat, celltype %in% c('Epithelial'))

# 赋一列
scRNA_chat$NMF_celltype='Epithelial_cells'
scRNA_chat=merge(scRNA_chat,scRNA_fibro.nmf)


meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))
#devtools::install_github("sqjin/CellChat")
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
                 weight.scale = T, label.edge= T, sources.use = 'Epithelial_cells',
                 title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, 
                 title.name = "Interaction weights/strength")


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = 'GAS', width = 8, height = 2.5, font.size = 10)


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2


gc()
### 转录因子
#####SCENIC#####
scRNAsub=readRDS('./scRNA_fibro_NMF.RDS')

# 有错记得本地安装
#devtools::install_github("aertslab/RcisTarget") 
#devtools::install_github("aertslab/SCENIC")
library(SCENIC)
library(Seurat)
#BiocManager::install(c("GENIE3","RcisTarget","CellChat"))

## 内存不够控制细胞数目！！！不要超过1000细胞！！！！！！！
set.seed(123456)
#a=sample(1:ncol(scRNAsub),200)
#scRNAsub@reductions$nmf=NULL
#scRNAsub=scRNAsub[,a]

exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo <-scRNAsub@meta.data

# 进入子文件夹
setwd("SCENIC") 

Idents(scRNAsub)=scRNAsub$NMF_celltype

### Initialize settings
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

Idents(scRNAsub) <- "NMF_cluster"

# 人物种
org='hgnc'


dbDir="./SCENIC" # RcisTarget databases location
myDatasetTitle="SCENIC Challenge" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
data(list=c("motifAnnotations_hgnc_v9","motifAnnotations_hgnc_v8"), package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org=org,
                                  dbDir=dbDir, 
                                  dbs = dbs)

scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"
??initializeScenic
# 节省内存，我们只取一个基因组文件
# 人：用下面！！！！！！！！！
scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr" = "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v9"
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
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123


exprMat_log <- log2(exprMat+1)
dim(exprMat)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# step2容易崩，节省内存只看top50
gc()
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50")) #** Only for toy run!!

gc()
# 回到单线程！否则容易报错！！！！！
scenicOptions@settings$nCores <- 1
library(foreach)

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="scenicOptions.Rds") # To save status

#============================================

#下次可以直接读取他,我们已经越过最耗时间的步骤
# 如果经常报错，不要重新readRDS和load，应不间断运行

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

fibro.scenic <- scRNAsub
# 结果导入endo.scenic
fibro.scenic@meta.data <- cbind(fibro.scenic@meta.data, aucell_regulonAUC.t[rownames(fibro.scenic@meta.data),])
dev.off()
DimPlot(fibro.scenic, reduction = "umap")
# 找6个
FeaturePlot(fibro.scenic, reduction = "umap", features = colnames(fibro.scenic@meta.data)[20:27], cols = c("yellow", "red"))

#Regulon scores heatmap 热图
Idents(fibro.scenic)=fibro.scenic$NMF_celltype
cells.ord.cluster <- fibro.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)

library(gplots)
library(pheatmap)

# 进一步标化
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(regulon.scores, 1, cal_z_score))
#pheatmap(data_subset_norm)

cluster.col <- data.frame(fibro.scenic@active.ident, row.names = names(fibro.scenic@active.ident))
colnames(cluster.col) <- "group"


#pheatmap::pheatmap(regulon.scores.scaled, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F, fontsize_row=5)
#pheatmap::pheatmap(data_subset_norm, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)


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


##### 关键CAF基因的表达
# 回到外面的文件夹!!!!!!!!!!!!!!!!!！否则必然报错
setwd('..')
gc()
scRNA=readRDS('./scRNA_fibro_NMF.RDS')

fibro_avg=AverageExpression(scRNA,group.by = 'NMF_celltype')
fibro_avg=as.data.frame(fibro_avg)
###  !!!!!!
colnames(fibro_avg)=c('C0','C1',"C2")
tme=read.table('TME_gene.txt')
tme=tme$V1

## ss必须把tme放前面！
ss=intersect(tme,rownames(fibro_avg))
rt=fibro_avg[ss,]

pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,cluster_rows = F,
                   color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100))


### CAF分型
### CAF亚型
library(data.table)
library(readxl)


gene_set = read.table(file='CAF_subtype.txt', header = T, sep = '\t',stringsAsFactors = F)
## 制作caf 的signature
list <- list()
for(i in 1:length(unique(gene_set$CAF_subtype))){
  list[[i]] <- gene_set$marker[gene_set$CAF_subtype== (unique(gene_set$CAF_subtype)[i])]
}
names(list)<- unique(gene_set$CAF_subtype)


scRNA=readRDS('./scRNA_fibro_NMF.RDS')
library(Seurat)
names(list)
scRNA=AddModuleScore(scRNA,features = list)

colnames(scRNA@meta.data)
## 注意选择列，你可能需要修改!!!!!!
df=scRNA@meta.data[,c(12:17)]


a=aggregate(df[,1:6],list(df$NMF_celltype),mean)
rownames(a)=a$Group.1
a$Group.1=NULL
a=a[,-1]
colnames(a)=names(list)
a=t(a)
a=a[,c(0,1)]
pheatmap::pheatmap(a,scale = 'row',cluster_cols = F,cluster_rows = F,
                   color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
)
