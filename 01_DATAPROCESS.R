rm(list=ls())
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

library(GEOquery)
library(Biobase)

setwd("01_数据预处理")

#读入ExpressionSet对象
gset51092=getGEO(filename = "GSE51092_series_matrix.txt",destdir = "./")
class(gset51092)
exprset=exprs(gset51092)
exprset[1:5,1:5]
#注释信息
GPL=getGEO(filename = 'GPL6884.soft',destdir = "./")
gpl=GPL@dataTable@table
ids=gpl[,c(1,6)]
ids<-ids[ids$ILMN_Gene != '', ]
rownames(ids)=ids[,1]
index<-ids[rownames(exprset),]
#保留最大表达量的探针对应的基因
index$max<-apply(exprset, 1, max)
identical(rownames(exprset),index$ID)
index<-index[order(index$ILMN_Gene,
                   index$max,
                   decreasing = T), ]
index<-index[!duplicated(index$ILMN_Gene), ]
exprset=exprset[index$ID,]
identical(rownames(exprset),index$ID)
rownames(exprset)=index$ILMN_Gene

#临床信息
phenoData<-pData(gset51092)
meta<-phenoData[, c("title","geo_accession")]

#制作分组信息
group<-strsplit(as.character(meta$title), " ")
group_list=sapply(group, "[",1)
meta$group<-group_list

write.csv(exprset,file ="exprset51092.csv", row.names = T, quote = F)
save(exprset,gpl,meta,file = 'gse51092.Rdata')

#表达矩阵和样本信息做好了，接下来检查一下
##1.判断是否需要normalization
##GEO的Data processing中显示已经过标准化：
##使用RMA对每个样本进行独立的标准化，之后进行log2转换和quantile标准化化，合并后使用ComBat处理批次效应
par(cex = 0.7)
rt<-exprset
n.sample=ncol(rt)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(rt, col = cols,main="expression value",las=2)

##2.判断是否需要log
expset<-rt
qx <- as.numeric(quantile(expset, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { expset[which(expset <= 0)] <- NaN
exprset<- log2(expset)
print("log2 transform finished")}else{print("log2 transform not needed")}
##不需要经过log2转换了###
save(exprset,meta,file = 'normalized_exprset51092.Rdata')
