#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("KEGG.db")
remotes::install_github("YuLab-SMU/createKEGGdb")
#创建自己的物种的包create_kegg_db，会自动创建名称为KEGG.db_1.0.tar,gz的包。物种名称的简写，在
library(createKEGGdb)
createKEGGdb::create_kegg_db("hsa")
#获取物种列表：'https://www.genome.jp/kegg/catalog/org_list.html'
#安装这个包(默认的包的路径在当前工作目录，根据实际情况修改路径)
install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")



#???ð?
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
#library(DO.db)
pvalueFilter=0.05      #pֵ????????
qvalueFilter=0.05      #????????pֵ????????

#devtools::install_github("YuLab-SMU/DOSE")
#devtools::install_github("YuLab-SMU/HDO.db")
#devtools::install_github('YuLab-SMU/clusterProfiler') 
library(DOSE)
library(HDO.db)
library(clusterProfiler)
#????ͼ?ε???ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")	
#setwd("C:\\Users\\lexb\\Desktop\\DRG\\28.KEGG")       #???ù???Ŀ¼
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #??ȡ?????ļ?

#??ȡ??????????????,??????????ת??Ϊ????id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #ȥ??????idΪNA?Ļ???
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

library(R.utils)

R.utils::setOption("clusterProfiler.download.method",'auto')
#KEGG????????
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1,use_internal_data = T)
?enrichKEGG
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter),]
#?????????????Ľ???
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#??????ʾͨ·????Ŀ
showNum=30     #??ʾ????????????ǰ30??ͨ·
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#??״ͼ
pdf(file="KEGG_barplot.pdf", width=8.5, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#????ͼ
pdf(file="KEGG_bubble.pdf", width=8.5, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

