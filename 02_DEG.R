rm(list=ls())
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

# install.packages("samr")

#数据输入
#表达矩阵exprSet
#样本分组group_list
library(limma)
library(samr)

setwd("../02_差异分析")

load("../01_数据预处理/normalized_exprset51092.Rdata")
exprSet<-exprset
group_list<-meta$group
#制作分组矩阵
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
#制作差异比较矩阵
contrast.matrix<-makeContrasts(contrasts=c('case-control'),levels = c("case","control"))
contrast.matrix

#开始差异分析
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
##default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
genes_list_limma=nrDEG[nrDEG$adj.P.Val < 0.05,]
deg_up=nrDEG[nrDEG$logFC>=1,]
deg_lo=nrDEG[nrDEG$logFC<=-1,]
table(duplicated(c(rownames(deg_up),rownames(deg_lo),rownames(genes_list_limma))))
save(genes_list_limma,deg_up,deg_lo,nrDEG,file = 'DEG_limma_exprset51092.Rdata')

topup<-deg_up[order(deg_up$adj.P.Val,
                    decreasing = F), ][1:21,]
toplo<-deg_lo[order(deg_lo$adj.P.Val,
                    decreasing = F), ]
topDEGs<-rbind(topup,toplo)
write.csv(topDEGs,file ="topDEGs_table1.csv", row.names = T, quote = F)

#SAM法
##samr函数
data=list(x=as.matrix(exprSet),y=as.numeric(as.factor(group_list)), 
          geneid=as.character(1:nrow(exprSet)),
          genenames=rownames(exprSet),
          logged2=TRUE)
samr.obj<-samr(data, resp.type="Two class unpaired", nperms=1000)#,fdr.output = 0.05

delta.table<-samr.compute.delta.table(samr.obj,nvals=100)
##计算不同delta(45度线到SAM定义阈值规则的上下平行线的垂直距离)与对应的fdr值
##这里筛选显著基因作为下一步WGCNA的input
##WGCNA本身不需要选择差异基因，这里主要考虑和作者的called数量保持一致
delta=1
##根据delta.table，考虑called基因的数量和fdr，选择SAM法的阈值delta
samr.plot(samr.obj,delta)
##可视化delta=1时筛选的显著基因
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
genes.up=siggenes.table$genes.up
genes.lo=siggenes.table$genes.lo
genes.lo_list=as.character(genes.lo[,2])
genes.up_list=as.character(genes.up[,2])
genes_list_samr=c(genes.lo_list,genes.up_list)
#d=unique(c(genes_list_samr,rownames(genes_list_limma)))
save(genes_list_samr,siggenes.table,samr.obj,delta.table,file = 'DEG_SAM_exprset51092.Rdata')
save(genes_list_samr,file = 'DEG_exprset51092.Rdata')
