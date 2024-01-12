# 数据处理
## 载入TCGA
load('./ESCA_tpm.Rdata')

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
gsva_matrix<- gsva(as.matrix(exprSet_tcga_mRNA), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)

data=gsva_matrix
#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

data=as.data.frame(t(data))

## 留下肿瘤
## 读入结果
tide_rt=read.csv('./TIDE_TCGA_output.csv',row.names = 1)
data=data[rownames(tide_rt),]
# 确认一致
tide_rt=tide_rt[rownames(data),]
identical(rownames(tide_rt),rownames(data))

data$Response=tide_rt$Responder

library(tidyverse)
rt=tidyr::pivot_longer(data,cols = -c('Response'),names_to = "Aggre_celltype",values_to = 'Abundance')


library(ggpubr)

ggpubr::ggboxplot(
  rt,
  x = "Aggre_celltype",
  y = "Abundance",
  color = "black",
  fill = "Response",
  xlab = "ICB Response",
  ylab = "Abundance", palette=c('#5d84a4','#dadada')
) +
  stat_compare_means(
    aes(group = Response),
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

# 逻辑回归
pFilter=1
rt=data
rt$Response=ifelse(rt$Response=='False',0,1)
outTab=data.frame()
for(gene in colnames(rt[,1:(ncol(rt)-1)])){
  set.seed(123456)
  glm=glm(Response ~ rt[,gene], data = rt,family= binomial(link='logit'))
  glmSummary = summary(glm)
  OR=exp(glm$coefficients)[2]
  OR_CI=exp(confint(glm,level = 0.95))
  OR.95L=OR_CI[2,1]
  OR.95H=OR_CI[2,2]
  glmP=glmSummary$coefficients[,"Pr(>|z|)"][2]
  if(glmP<pFilter){
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       OR=OR,
                       OR.95L=OR.95L,
                       OR.95H=OR.95H,
                       pvalue=glmP))
  }
}

write.table(outTab,file="logistic_TIDE_tcga.txt",sep="\t",row.names=F,quote=F)


####acrg
# 数据处理
## 载入TCGA
load('GSE53625.Rdata')
#TIDE <- exprSet#行是基因，列是样本
library(data.table)
pdata = fread("./GSE53624suv.txt",data.table = F,sep = '\t')
rownames(pdata)=pdata[,1]
pdata=pdata[,-1]
#TCGA需要挑样本，GEO全是肿瘤则不需要
ss=intersect(rownames(pdata),colnames(exprSet))
exprSet=exprSet[,ss]

hc=grep('normal',pdata$title)
length(hc)
pdata=pdata[-hc,]
length(rownames(pdata))
exprSet=exprSet[,-hc]

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

data=gsva_matrix


data=as.data.frame(t(data))

## 留下肿瘤
## 读入结果
tide_rt=read.csv('./TIDE_GEO_output.csv',row.names = 1)
data=data[rownames(tide_rt),]
# 确认一致
tide_rt=tide_rt[rownames(data),]
identical(rownames(tide_rt),rownames(data))

data$Response=tide_rt$Responder

library(tidyverse)
rt=tidyr::pivot_longer(data,cols = -c('Response'),names_to = "Aggre_celltype",values_to = 'Abundance')


library(ggpubr)

ggpubr::ggboxplot(
  rt,
  x = "Aggre_celltype",
  y = "Abundance",
  color = "black",
  fill = "Response",
  xlab = "ICB Response",
  ylab = "Abundance", palette=c('#5d84a4','#dadada')
) +
  stat_compare_means(
    aes(group = Response),
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

# 逻辑回归
pFilter=1
rt=data
rt$Response=ifelse(rt$Response=='False',0,1)
outTab=data.frame()
for(gene in colnames(rt[,1:(ncol(rt)-1)])){
  set.seed(123456)
  glm=glm(Response ~ rt[,gene], data = rt,family= binomial(link='logit'))
  glmSummary = summary(glm)
  OR=exp(glm$coefficients)[2]
  OR_CI=exp(confint(glm,level = 0.95))
  OR.95L=OR_CI[2,1]
  OR.95H=OR_CI[2,2]
  glmP=glmSummary$coefficients[,"Pr(>|z|)"][2]
  if(glmP<pFilter){
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       OR=OR,
                       OR.95L=OR.95L,
                       OR.95H=OR.95H,
                       pvalue=glmP))
  }
}

write.table(outTab,file="logistic_TIDE_acrg.txt",sep="\t",row.names=F,quote=F)


#### 组合OR值的气泡热图
#可以组合更多数据集
rt1=read.table('./logistic_TIDE_tcga.txt',header = T)
rt2=read.table('./logistic_TIDE_acrg.txt',header = T)
rt1$Dataset='TCGA'
rt2$Dataset='GEO'


rt=rbind(rt1,rt2)
rt$`-logP`=-log(rt$pvalue)
rt$logOR=log(rt$OR)
library(ggplot2)
rt$logOR=ifelse(rt$logOR>3,3,rt$logOR)
rt$logOR=ifelse(rt$logOR< -3,-3,rt$logOR)

ggplot(data=rt)+
  geom_point(aes(y=Dataset,x=gene,fill=logOR,size= `-logP`),
             color='black',shape=21,stroke=1.5)+
  scale_fill_gradientn(colours = c('#403D76','#E3B635','#C02E20'),limits=c(-3,3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,hjust=1),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  scale_size_area(breaks=c(1,2,3))


load('./expcli_IMvigor210.Rdata')
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
gsva_matrix<- gsva(as.matrix(eset), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)

data=gsva_matrix
#删掉正常样品
data=as.data.frame(t(data))

identical(rownames(pdata),rownames(data))

data$Response=pdata$binaryResponse
data=data[,-13]
library(tidyverse)
rt=tidyr::pivot_longer(data,cols = -c('Response'),names_to = "Aggre_celltype",values_to = 'Abundance')

rt=na.omit(rt)

rt$Response=factor(rt$Response,levels = c('SD/PD','CR/PR'))
library(ggpubr)

ggpubr::ggboxplot(
  rt,
  x = "Aggre_celltype",
  y = "Abundance",
  color = "black",
  fill = "Response",
  xlab = "ICB Response",
  ylab = "Abundance", palette=c('#5d84a4','#dadada')
) +
  stat_compare_means(
    aes(group = Response),
    label = "p.signif", 
    #method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))

data$futime=pdata$os
data$fustat=pdata$censOS

colnames(rt)
### 反复重命名是为了少修改代码
rt=data
colnames(data)
library(survival)
library(survminer)
res.cut=surv_cutpoint(rt, time="futime", 
                      event="fustat",
                      variables='Non-Aggre-B_cells-C2')
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,'Non-Aggre-B_cells-C2']<= cutoff, "Low", "High")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c('#ebc03e','#b42e20')
bioCol=bioCol[1:length]
ggsurvplot(fit, 
           data=data,
           conf.int=F,
           pval=pValue,
           pval.size=6,
           legend.title='Non-Aggre-B_cells-C2',
           legend.labs=levels(factor(data[,"group"])),
           legend = c(0.7, 0.8),
           font.legend=12,
           xlab="Time(Month)",
           palette = bioCol,
           #surv.median.line = "hv",
           risk.table=F,
           cumevents=F,
           risk.table.height=.25)

