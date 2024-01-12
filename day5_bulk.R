# 数据处理
## 载入TCGA
load('./ESCA_tpm.Rdata')

## Bulk层面的聚集体自噬评估
gene=read.table('REACTOME_AGGREPHAGY.v7.5.1.gmt')
gene=gene[,3:ncol(gene)]
gene=t(gene)

gene_set=list(gene)
names(gene_set)='AGGREPHAGY'


library(genefilter)
library(GSVA)
library(Biobase)
library(edgeR)
# gsva方法
gsva_matrix<- gsva(as.matrix(exprSet_tcga_mRNA), gene_set,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)

## 肿瘤和正常
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))

rt=as.data.frame(t(gsva_matrix))
rt$group=group_list


library(ggpubr)

library(tidyverse)
colnames(rt)

ggboxplot(rt,x = "group",
          y = "AGGREPHAGY",
          color = "black",
          fill = "group",
          xlab = "group",
          ylab = "AGGREPHAGY", palette=c('#b42e20','#ebc03e')
)+stat_compare_means()


## 整理好之前获取的矩阵
gene_set = read.table(file='aggre_celltype.txt', header = T, sep = '\t',stringsAsFactors = F)
#names(gene_set)<-c("Aggre_celltype","marker")
list <- list()
for(i in 1:length(unique(gene_set$Aggre_celltype))){
  list[[i]] <- gene_set$marker[gene_set$Aggre_celltype== (unique(gene_set$Aggre_celltype)[i])]
}
names(list)<- unique(gene_set$Aggre_celltype)

library(genefilter)
library(GSVA)
library(Biobase)
library(edgeR)
gsva_matrix<- gsva(as.matrix(exprSet_tcga_mRNA), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)

## 肿瘤和正常
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))

rt=as.data.frame(t(gsva_matrix))
rt$group=group_list

library(ggpubr)

library(tidyverse)
# 变长矩阵
rt=tidyr::pivot_longer(rt,cols = -c('group'),names_to = "Aggre_celltype",values_to = 'Abundance')


library(ggpubr)

ggboxplot(
  rt,
  x = "Aggre_celltype",
  y = "Abundance",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "Abundance", palette=c('#b42e20','#ebc03e')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))



### 生存分析
data=gsva_matrix
#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

data=as.data.frame(t(data))

# 最优Cutoff寻找预后影响的新型聚集体自噬免疫细胞
## 读取生存数据
suv=read.table('./TCGA-ESCA.survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")

## 数据合并
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
#你需要牢记这样的结构，用于cox
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
## 写出聚集体自噬免疫细胞的基础矩阵
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)


## 单因素coX
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件

library(survminer)
library(survival)
library(dplyr)
library(caret)

#单因素COX
pFilter=1                                                #显著性过滤标准     #设置工作目
outTab=data.frame()
sigGenes=c("futime","fustat")

#循环
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="uniCox_tcga.txt",sep="\t",row.names=F,quote=F)


rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件


## K-M
rt$futime=rt$futime/30
immune_p=c()
immune_figure=list()
library(survival)
library(survminer)

#install.packages('cowplot')
library(cowplot)
dev.off()
dir.create('k_m')
for (i in rownames(gsva_matrix)) {
  res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(data[,i]<= cutoff, "Low", "High")
  data=rt
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = data)
  length=length(levels(factor(data[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  immune_p=c(immune_p,pValue)
  if(pValue<0.05){
    fit <- survfit(Surv(futime, fustat) ~ group, data = data)
    bioCol=c('#ebc03e','#b42e20')
    bioCol=bioCol[1:length]
    p=ggsurvplot(fit, 
                 data=data,
                 conf.int=F,
                 pval=pValue,
                 pval.size=6,
                 legend.title=i,
                 legend.labs=levels(factor(data[,"group"])),
                 legend = c(0.7, 0.8),
                 font.legend=12,
                 xlab="Time(Months)",
                 palette = bioCol,
                 #surv.median.line = "hv",
                 risk.table=F,
                 cumevents=F,
                 risk.table.height=.25)
    ggsave2(filename = paste0('./k_m/',i,'.pdf'),width = 4,height = 4)
  }
}


#######################ACRG：另一个胃癌数据库################

load('./GSE53625.Rdata')
exprSet_acrg=exprSet
gene_set = read.table(file='aggre_celltype.txt', header = T, sep = '\t',stringsAsFactors = F)
gene_set=gene_set[,-3]
colnames(gene_set)=c("Aggre_celltype","marker")
list <- list()
for(i in 1:length(unique(gene_set$Aggre_celltype))){
  list[[i]] <- gene_set$marker[gene_set$Aggre_celltype== (unique(gene_set$Aggre_celltype)[i])]
}
names(list)<- unique(gene_set$Aggre_celltype)

library(genefilter)
library(GSVA)
library(Biobase)
library(edgeR)

gsva_matrix<- gsva(as.matrix(exprSet_acrg), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)

### 另一个胃癌数据集生存分析
data=gsva_matrix
data=as.data.frame(t(data))

# 最优Cutoff寻找预后影响的新型免疫细胞
## 读取生存数据
suv=read.table('./GSE53624suv.txt',row.names = 1,header = T,check.names = F)


cli=dplyr::select(suv,'futime','fustat')
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/30
## 数据合并并输出结???
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
#你需要牢记这样的结构，用于cox
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
## 写出GPR的基础矩阵
write.table(out,file="expTime_acrg.txt",sep="\t",row.names=F,quote=F)





## 单因素coX------------------------
rt=read.table("expTime_acrg.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件

library(survminer)
library(survival)
library(dplyr)
library(caret)

#单因素COX
pFilter=1                                                #显著性过滤标准     #设置工作目
outTab=data.frame()
sigGenes=c("futime","fustat")

#循环
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="uniCox_acrg.txt",sep="\t",row.names=F,quote=F)


#########------------------
rt=read.table("expTime_acrg.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件

## K-M
immune_p=c()
immune_figure=list()
library(survival)
library(survminer)

#install.packages('cowplot')
library(cowplot)
dir.create('k_m_acrg')
for (i in rownames(gsva_matrix)) {
  res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(data[,i]<= cutoff, "Low", "High")
  data=rt
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = data)
  length=length(levels(factor(data[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  immune_p=c(immune_p,pValue)
  if(pValue<0.05){
    fit <- survfit(Surv(futime, fustat) ~ group, data = data)
    bioCol=c('#ebc03e','#b42e20')
    bioCol=bioCol[1:length]
    p=ggsurvplot(fit, 
                 data=data,
                 conf.int=F,
                 pval=pValue,
                 pval.size=6,
                 legend.title=i,
                 legend.labs=levels(factor(data[,"group"])),
                 legend = c(0.7, 0.8),
                 font.legend=12,
                 xlab="Time(Months)",
                 palette = bioCol,
                 #surv.median.line = "hv",
                 risk.table=F,
                 cumevents=F,
                 risk.table.height=.25)
    ggsave2(filename = paste0('./k_m_acrg/',i,'.pdf'),width = 4,height = 4)
  }
}



#### 组合HR值的气泡热图
#可以组合更多数据集
rt1=read.table('./uniCox_tcga.txt',header = T)
rt2=read.table('./uniCox_acrg.txt',header = T)
rt1$Dataset='TCGA'
rt2$Dataset='GEO'


rt=rbind(rt1,rt2)
rt$`-logP`=-log(rt$pvalue)
rt$logHR=log(rt$HR)
library(ggplot2)
rt$logHR=ifelse(rt$logHR>3,3,rt$logHR)
rt$logHR=ifelse(rt$logHR< -3,-3,rt$logHR)

ggplot(data=rt)+
  geom_point(aes(y=Dataset,x=id,fill=logHR,size= `-logP`),
             color='black',shape=21,stroke=1.5)+
  scale_fill_gradientn(colours = c('#403D76','#E3B635','#C02E20'),limits=c(-3,3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,hjust=1),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  scale_size_area(breaks=c(1,2,3))
  


