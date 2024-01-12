usethis::browse_github_pat()
usethis::edit_r_environ()
#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
devtools::install_github("YuLab-SMU/DOSE")
devtools::install_github("YuLab-SMU/HDO.db")
devtools::install_github('YuLab-SMU/clusterProfiler') #作者：生信私学 https://www.bilibili.com/read/cv22520449/ 出处：bilibili
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
for (i in c('TUBA1A+B_cell-C0','UBE2N+B_cell-C1','Non-Aggre-B_cells-C2')) {
  a <- gene %>% 
    filter(Aggre_celltype == i)
  genelist <- bitr(a$marker,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  gene_f <- a %>% 
    right_join(genelist, by = c("marker"="SYMBOL" ))
  kegg[[i+1]] <- enrichKEGG(gene =  gene_f$ENTREZID,
                            organism     = 'hsa',
                            pvalueCutoff = 0.1,
                            qvalueCutoff =0.1)
}

#
####!!!!
names(kegg) <- paste0('NMF_cluster',c(0,1,2,3,4,5))

## 如果有NULL的无富集结果的情况，需要去除！！！
# 例如
kegg$NMF_cluster4=NULL
### 有内容的话不去除

celltype_pathway <- data.frame()

##### ！！！！！！！！！！！！！！！！！
for (i in 1:5) {
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
celltype_pathway <- celltype_pathway[c(1:3,46:48,74,75:77,87),]

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




