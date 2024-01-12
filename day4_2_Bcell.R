####################################B细胞##########################################

gc()
load('./scRNA_anno.rds')
library(Seurat)
DimPlot(scRNA,split.by = 'orig.ident',reduction = 'tsne')

# B胞的NMF和分析
scRNA_Bcell=subset(scRNA,celltype %in% 'B_cells')

dim(scRNA_Bcell)

## 重新聚类分群
scRNA_Bcell<- ScaleData(scRNA_Bcell)

scRNA_Bcell<- RunPCA(scRNA_Bcell) 
plot1 <- DimPlot(scRNA_Bcell, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA_Bcell, ndims=20, reduction="pca") 
plotc <- plot1+plot2
plotc

#PC选取比较随意
pc.num=1:15

scRNA_Bcell <- FindNeighbors(scRNA_Bcell, dims = pc.num) 

scRNA_Bcell <- FindClusters(scRNA_Bcell)
table(scRNA_Bcell@meta.data$seurat_clusters)

#tSNE可视化
scRNA_myelo = RunTSNE(scRNA_Bcell, dims = pc.num)
plot1 <-DimPlot(scRNA_Bcell, reduction = "tsne",label = T,split.by = 'orig.ident') 

#UMAP可视化
scRNA_myelo <- RunUMAP(scRNA_Bcell, dims = pc.num)
plot2 <-DimPlot(scRNA_Bcell, reduction = "umap",label=T,split.by = 'orig.ident') 
plotc <- plot1+plot2
plotc
habermann_imm <- c("JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", 'CD20',"CD79A")
## 可以看到1，2是巨噬细胞，0是肥大细胞，3是单核细胞
library(ggplot2)
DotPlot(scRNA_Bcell, features = habermann_imm,group.by = "seurat_clusters") +
  coord_flip()+scale_colour_gradientn(colours =c('#dadada','#bc3c29'))

## 
B_subtype <- c('B_cells','B_cells','B_cells','Plasma_cells','B_cells','Plasma_cells','Plasma_cells','B_cells')

Idents(scRNA_Bcell) <- scRNA_Bcell@meta.data$seurat_clusters
names(B_subtype) <- levels(scRNA_Bcell)
scRNA_Bcell<- RenameIdents(scRNA_Bcell, B_subtype)

scRNA_Bcell@meta.data$B_subtype <- Idents(scRNA_Bcell)

Idents(scRNA_Bcell)=scRNA_Bcell@meta.data$B_subtype

DimPlot(scRNA_Bcell, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                      '#7bc7cd','#5d84a4','#dadada'),
        reduction = 'tsne',split.by = 'orig.ident') 



## 取出肿瘤
scRNA_Bcell=subset(scRNA_Bcell,orig.ident %in% 'Tumor')


gene=read.table('REACTOME_AGGREPHAGY.v7.5.1.gmt')
gene=gene[,3:46]
gene=t(gene)

# 聚集体自噬aa
scRNA_Bcell$id=colnames(scRNA_Bcell)
scRNA_Bcell_aa=scRNA_Bcell[rownames(scRNA_Bcell) %in% gene,]

df <- scRNA_Bcell_aa@assays$RNA@data
df=as.data.frame(df)

df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]

## 一点没表达的细胞去除
scRNA_Bcell=subset(scRNA_Bcell,id %in% colnames(df))

# nmf
library(NMF)
res <- nmf(df, 10, method = "snmf/r", seed = 'nndsvd')

## NMF结果返回suerat对象
scRNA_Bcell@reductions$nmf <- scRNA_Bcell@reductions$pca
scRNA_Bcell@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_Bcell@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类!重要
set.seed(999)

scRNA_Bcell.nmf <- RunUMAP(scRNA_Bcell, reduction = 'nmf', dims = 1:10) 
scRNA_Bcell.nmf <- FindNeighbors(scRNA_Bcell.nmf,reduction = 'nmf', dims = 1:10)
scRNA_Bcell.nmf <- FindClusters(scRNA_Bcell.nmf)
scRNA_Bcell.nmf$NMF_cluster=scRNA_Bcell.nmf$seurat_clusters


## 结果可视化  
DimPlot(scRNA_Bcell.nmf, label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                        '#7bc7cd','#5d84a4','#dadada')) 


df=FindAllMarkers(scRNA_Bcell.nmf,logfc.threshold = 0.5,only.pos = T)


write.csv(df,file ='deg_B_cell.csv',quote=F)


library(tidyverse)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
library(R.utils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
gene <- read.csv("deg_B_cell.csv")
kegg <- list()

#修改!!!!!
for (i in 0:5) {
  a <- gene %>% 
    filter(cluster == i)
  genelist <- bitr(a$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  gene_f <- a %>% 
    right_join(genelist, by = c("gene"="SYMBOL" ))
  kegg[[i+1]] <- enrichKEGG(gene =  gene_f$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.1,
                            qvalueCutoff =0.1)
}

#
####!!!!
names(kegg) <- paste0('NMF_cluster',0:5)

## 如果有NULL的无富集结果的情况，需要去除！！！
# 例如
kegg$NMF_cluster4=NULL
### 有内容的话不去除

celltype_pathway <- data.frame()

##### ！！！！！！！！！！！！！！！！！
for (i in 1:6) {
  kegg_single <- kegg[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column("name")
  kegg_single<-  kegg_single %>% 
    mutate(log.p.adjust = -log( kegg_single$p.adjust)) %>% 
    mutate(celltype_pathway = names(kegg)[i])
  celltype_pathway <- rbind(celltype_pathway, kegg_single)    #合并                 
}

celltype_pathway=celltype_pathway[celltype_pathway$p.adjust<0.05,]

## 去除行名，使其真的为数字编号的行名
rownames(celltype_pathway)=NULL
celltype_pathway <-celltype_pathway %>% 
  dplyr::select(5,13,14) 

## 每个cluster找前3条
celltype_pathway <- celltype_pathway[c(1:3,29:31,34:36,52:54,73:75,85:87),]

library(RColorBrewer)
library(ggplot2)
celltype_pathway$Description <- factor(celltype_pathway$Description)

ap <- ggplot(celltype_pathway,aes(x=celltype_pathway,y=Description))+
  geom_tile(aes(fill=log.p.adjust),color = "grey")+
  scale_fill_gradient(low = "#F7ED7A",high = "#E03220")+
  theme_classic()+theme(axis.text.x = element_text(angle=90))+
  labs(title = "KEGG Pathways",
       x = element_blank(),
       y=element_blank()
  )

ap





# test some genes
FeaturePlot(scRNA_Bcell.nmf,features = 'TUBB4B',label = T)
DoHeatmap(scRNA_Bcell.nmf,features = rownames(df),assay = 'RNA',slot = 'count')

#'Non-Aggre-B_cells-C1'
## 定义新的细胞亚群！注意修改代码！！！！！！！！！！！！！！
NMF_celltype <- c('TUBA1A+B_cell-C0', 'UBE2N+B_cell-C1',"TUBA1A+B_cell-C0","TUBA1A+B_cell-C0","Non-Aggre-B_cells-C2",'Non-Aggre-B_cells-C2')

Idents(scRNA_Bcell.nmf) <- scRNA_Bcell.nmf@meta.data$NMF_cluster
names(NMF_celltype) <- levels(scRNA_Bcell.nmf)
scRNA_Bcell.nmf<- RenameIdents(scRNA_Bcell.nmf, NMF_celltype)

scRNA_Bcell.nmf@meta.data$NMF_celltype <- Idents(scRNA_Bcell.nmf)

Idents(scRNA_Bcell.nmf)=scRNA_Bcell.nmf@meta.data$NMF_celltype

DimPlot(scRNA_Bcell.nmf,group.by = 'NMF_celltype',
        cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
               '#7bc7cd','#5d84a4'),
        label = T,label.size = 3)

saveRDS(scRNA_Bcell.nmf,file ='scRNA_Bcell_NMF.RDS')

proportion=scRNA_Bcell.nmf@meta.data
colnames(proportion)
proportion=proportion[,c('B_subtype','NMF_celltype')]

library(ggstatsplot)
ggstatsplot::ggbarstats(proportion,x='NMF_celltype',y='B_subtype')


scRNA_Bcell.nmf<-scRNA_Bcell_NMF
########## 细胞通讯(B细胞和T细胞）#################################
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')
scRNA_chat <-subset(scRNA_chat, celltype=='T_cells')

# 赋一列
scRNA_chat$NMF_celltype='T_cells'
scRNA_chat=merge(scRNA_chat,scRNA_Bcell.nmf)

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
                 weight.scale = T, label.edge= T,targets.use = 'T_cells',
                 title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, 
                 title.name = "Interaction weights/strength")

##### 通讯集合
scRNA_fibro=readRDS('scRNA_fibro_NMF.RDS')
scRNA_T=readRDS('scRNA_T_NMF.RDS')
scRNA_B=readRDS('scRNA_Bcell_NMF.RDS')
scRNA_mac=readRDS('scRNA_mac_NMF.RDS')


########## 细胞通讯(B细胞和上皮）#################################
scRNA=readRDS('scRNA_anno.rds')
scRNA_chat <-subset(scRNA, orig.ident =='Tumor')
scRNA_chat <-subset(scRNA_chat, celltype=='Epithelial')

# 赋一列
scRNA_chat$NMF_celltype='Epithelial_cells'
scRNA_chat=merge(scRNA_chat,c(scRNA_fibro,
                              scRNA_T,scRNA_B,scRNA_mac))



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
                 weight.scale = T, label.edge= T,
                 title.name = "Number of interactions")
dev.off()

library(viridis)

rt=df.net
colnames(rt)
p1 <- rt %>% 
  ggplot(aes(x = target,y = interaction_name,
             fill = prob, 
             split="source"),
         cex.axis=0.5) +
  geom_tile() +
  theme_bw() +
  # Because we need the x and y axis to display every node,
  # not just the nodes that have connections to each other,
  # make sure that ggplot does not drop unused factor levels
  #  scale_x_discrete(drop = FALSE) +
  #  scale_y_discrete(drop = FALSE) +
  theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0,vjust = 0.8),
    # Force the plot into a square aspect ratio
    aspect.ratio = 6
    # 去除图例
    #legend.position = "none"
  )
pdf("imm_hot_clauster.pdf",width =30,height = 12)
print(p1+labs(y="Interaction")+scale_fill_viridis()+facet_grid(.~source)+xlab(NULL))
dev.off()
