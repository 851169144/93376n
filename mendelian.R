## 关键基因的孟德尔随机化分析
library(dplyr)
# 读取数据
scRNA=readRDS('./scRNA_anno.RDS')
scRNAsub=readRDS('./scRNA_T.RDS')
library(Seurat)

scRNA_other=subset(scRNA,celltype !='T_cells')
scRNA_CM=subset(scRNAsub,T_celltype == 'CD8_CM')

##！！ 加一列
gc()
scRNA_other$T_celltype=scRNA_other$celltype

scRNA_compare=merge(scRNA_other,scRNA_CM)
table(scRNA_compare$T_celltype)

gc()
rm()


## 1. CM相对其他T亚群
df_CM=FindMarkers(scRNAsub,ident.1 = 'CD8_CM',only.pos = T,logfc.threshold = 0.5)

## 2. CM相对其他细胞
df_T=FindMarkers(scRNA_compare,ident.1 = 'CD8_CM',only.pos = T,logfc.threshold = 0.5)

ss=intersect(rownames(df_CM),rownames(df_T))

ss

## 保存，后面要用
save(ss,file ='key_marker_gene.Rdata')

load('key_marker_gene.Rdata')
write.csv(ss,file ='key_marker_gene.csv',quote = F)

## 自助教程：
## marker基因的富集分析怎么做,负值csv中的基因列到网站，参见下面教程
#https://mp.weixin.qq.com/s/ClHOFvw3GSM9wvmIPip4VA

gene=ss

gene=read.csv(file ='aggre_celltype.csv',header = T)
gene=gene$marker[!duplicated(gene$marker)]
## 转ENSEMBL
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler",force = TRUE)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
df=bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
id=df$ENSEMBL
id=id[!duplicated(id)]
id



## 载入包
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#https://gwas.mrcieu.ac.uk/datasets/

# 暴露数据的SNP，工具变量
exposure_id=paste0('eqtl-a-',id)
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)

R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)

exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)

exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]

exposure_dat$exposure=stringr::str_sub(exposure_dat$exposure,1,15)

exposure_dat$ENSEMBL =exposure_dat$exposure

library(tidyverse)
exposure_dat= left_join(exposure_dat,df, by = "ENSEMBL")

exposure_dat$exposure=NULL
exposure_dat$exposure=exposure_dat$SYMBOL
exposure_dat$SYMBOL=NULL
save(exposure_dat,file ='exposure_dat.Rdata')

load('exposure_dat.Rdata')



rm(scRNA_CM)
rm(scRNA_compare)
rm(scRNA)
gc()


#提取结局数据#ieu-b-4960 #bbj-a-117
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="bbj-a-117")

# 取交集
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]


harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)



## MR分析

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}




mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)

save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='mr_input_res.Rdata')

load('mr_input_res.Rdata')

#  异质性 ，用处不大
# 当存在显著的异质性（p<0.05）时，结果不太可靠
heter_dat=mr_heterogeneity(harmonised_dat)
write.csv(heter_dat, file="table.heterogeneity.csv", row.names=F)

# 水平多效性检验，用处不大
#如果该截距项与0差异很大（pval < 0.05），说明存在水平多效性
pleio_dat=mr_pleiotropy_test(harmonised_dat)
write.csv(pleio_dat, file="table.pleiotropy.csv", row.names=F)


### table
table1 <- mr_res %>% 
  filter(pval < 0.05,
         method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  left_join(exposure_dat, by = "exposure")

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = exposure, `ENSEMBL ID` = ENSEMBL,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                PVE, `F statistics`)

## 保存，后面常用
save(table1,file ='table1.Rdata')



### 可视化
volcano_plot <- function(.data, 
                         number_comparasion = 1,
                         title = "eQTL",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         legend.position = "none")
{
  p_thershold <- 0.05/number_comparasion
  
  p <- .data %>% 
    rename(beta := !!col_beta,
           size := !!col_size,
           label := !!col_label) %>% 
    mutate(x = beta,
           y = -log10(pval),
           label = ifelse(pval < p_thershold, label, NA)) %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = size), alpha = 0.5, color = "#0072b5") +
    geom_vline(xintercept = 0, linetype = 2)+
    geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.position = legend.position)+
    labs(x = "ln(OR)", 
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title) +
    scale_size(name = "PVE",
               breaks = c(0.2*1:3)) +
    ggrepel::geom_label_repel(aes(label = label),size = 3)
  plot(p)
}



  mr_res %>% 
    filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
    volcano_plot(number_comparasion = 1)



## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl！！
# 取一个有多个工具变量的基因
mr_res1=mr_res[mr_res$exposure=='VSIG4',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='VSIG4',]

#绘制散点图
#每一个点其实代表的就是一个IV，每个点上的线实际反映的是95%置信区间，横坐标是SNP对暴露的效应，纵坐标是SNP对结局的效应
mr_scatter_plot(mr_res1, harmonised_dat1)

#森林图 每一条水平实线反映的是单个SNP利用Wald ratio方法估计出来的结果
res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

#漏斗图  这和我们之前发现这些SNP间的异质性很大的结果相符。
mr_funnel_plot(singlesnp_results = res_single)

#留一法敏感性分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))


mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"OR.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
      )


