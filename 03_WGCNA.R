rm(list=ls())
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

setwd("../03_WGCNA")

load("../01_数据预处理/normalized_exprset51092.Rdata")
load("../02_差异分析/DEG_exprset51092.Rdata")
#WGCNA
exprset<-exprset[genes_list_samr,]
datExpr<-as.data.frame(t(exprset))
datTraits<-meta
identical(rownames(datExpr),rownames(datTraits))
save(datExpr,datTraits,file ="./WGCNA_input.Rdata")

library(PerformanceAnalytics)
library(WGCNA)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

rm(list=ls())
Sys.setenv(LANGUAGE = "en") 
load("./WGCNA_input.Rdata")

##相关性计算，构建无标度网络需要考虑的几个参数设置：
##type = "unsigned" #官方推荐 "signed" 或 "signed hybrid"
corType = "pearson" #官方推荐 biweight mid-correlation & bicor
##corFnc = ifelse(corType=="pearson", cor, bicor)
##maxPOutliers = ifelse(corType=="pearson",1,0.05) #选择"bicor"需要设置（对二元变量，如样本性状信息计算相关性时，或基因表达严重依赖于疾病状态时）
##robustY = ifelse(corType=="pearson",T,F) #关联样品性状的二元变量时设置

# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的R方和平均连接度，RsquaredCut为期望的R^2阈值
sft = pickSoftThreshold(datExpr,powerVector = powers,RsquaredCut = 0.9,verbose = 5)
# 作图：
pdf(file = "Fig1A.pdf",width = 9,height = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
softPower<-sft$powerEstimate
softPower
#这里计算的阈值为7，与原文的8不一致，下面以7继续进行分析
#可能由于所选统计学方法不同或差异分析最终筛选到的基因有差异

# 一步法构建共表达矩阵
# 使用上一步获得的sft$powerEstimate=7,设置mergeCutHeight=0.2
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned", 
  minModuleSize = 10,#设置每个基因模块最少的基因数目为10。
  reassignThreshold = 0, mergeCutHeight = 0.2,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)
## saving TOM for block 1 into file blockwiseTOM-block.1.RData
table(net$colors)

cor<-stats::cor
#层级聚类树展示模块 
#这里用不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任何模块的那些基因，
#如果灰色模块里面的基因太多，那么前期对表达矩阵挑选基因的步骤可能就不太合适。
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf(file = "dendrogram.pdf",width = 6,height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.
dev.off()

##动态剪接法的聚类树和合并后的聚类树
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
# 计算基因之间的相异度
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = 10
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
mergedColors = labels2colors(net$colors)
table(mergedColors)
pdf("Fig1B.pdf",width = 8,height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##可视化基因网络（TOM plot）
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^softPower
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
pdf("Fig1C_network_heatmap.pdf",width = 9,height = 9)
# 这一部分比较耗时，行和列同时做层级聚类
TOMplot(plotTOM, geneTree, mergedColors, col=colorRampPalette(colors = c("red","yellow","#FFFACD"))(50),
        main = "Network heatmap plot, all genes")
dev.off()

# Calculate eigengenes，官方是重新计算的，这里我们改下列名字就好
# MEList = moduleEigengenes(datExpr, colors = dynamicColors)
# MEs = MEList$eigengenes
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs = orderMEs(MEs_col)
### 模块与表型数据关联
design=model.matrix(~0+ datTraits$group)
design=as.matrix(design)
colnames(design)=c("SS","control")
moduleColors <- labels2colors(net$colors)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("Fig2A_Module-trait_relationships.pdf",width = 6,height = 8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               xLabelsAngle = 0,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               #showCols = 1
)
##原文只展示了"SS",在WGCNAv1.68及以前的版本中可以实现，1.69之后的版本不允许展示1个Col了
dev.off()
table(moduleColors)

### 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

#性状与每个基因表达量相关性在各个模块的均值作为该性状在该模块的显著性，显著性最大模块与该性状最相关
GS1 <- as.numeric(WGCNA::cor(design[,1],datExpr,use="p",method="pearson"))
# 显著性绝对值：
GeneSignificance <- abs(GS1)
# 获得该性状在每个模块中的显著性：
ModuleSignificance <- tapply(GeneSignificance,moduleColors,mean,na.rm=T)
ModuleSignificance <- as.matrix(t(ModuleSignificance))
sd <- tapply(GeneSignificance,moduleColors,sd,na.rm=T)
SE.mean <- sd/sqrt(table(moduleColors))
pdf("Fig2B_module_signif.pdf",width = 8,height = 6)
barplot_SS <- barplot(ModuleSignificance[1,],names.arg=F,
        ylab = "Gene Significance",ylim = c(0,0.5),
        col = colnames(ModuleSignificance),)
arrows(barplot_SS,ModuleSignificance[1,]+SE.mean,barplot_SS,ModuleSignificance[1,]-SE.mean,length =0.05, angle = 90, code = 3)
dev.off()
# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr, design, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr, design, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "turquoise"
pheno = "SS"
modNames = substring(colnames(MEs), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module
pdf("Fig2C_SS_turquoise.pdf",width = 8,height = 8)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##选择hub基因
##2个标准：gene significance >0.4, Module membership value > 0.9
##实际操作结合原文的fig2c，这里选择0.35&0.8
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
# 找到满足条件的基因：
hubGene <- (GeneSignificance > 0.35&(abs(datKME[,paste("MM.",module,sep="")]))>0.8)
hub=colnames(datExpr)
hub=hub[hubGene==TRUE&moduleGenes==TRUE]
hub
save(hub,file='hubgene.Rdata')

# hub 基因热图：
pdf("Fig2D_signed_correlations.pdf",width = 8,height = 8)
plotNetworkHeatmap(datExpr,
                   plotGenes = hub,
                   networkType = "signed",
                   useTOM = TRUE,
                   power=softPower,
                   main="signed correlations")
dev.off()

##对moduleGenes进行功能注释
##原文用的DAVID做GO富集，GOplot绘图，这里直接用clusterProfiler了
moduleGenes=colnames(datExpr)[moduleGenes==TRUE]
save(moduleGenes,file='GOanalysis.Rdata')

tmp <- bitr(moduleGenes, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
moduleGenes <- tmp$ENTREZID
egoBP <- enrichGO(gene = moduleGenes,  #基因列表文件中的基因名称
                  OrgDb = org.Hs.eg.db,  #指定物种的基因数据库
                  keyType = 'ENTREZID',  #指定给定的基因名称类型
                  ont = 'BP',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                  pAdjustMethod = 'BH',  #指定 p 值校正方法
                  pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                  readable = TRUE)
egoBP1 <- clusterProfiler::simplify(egoBP,measure = "Wang",
                                    cutoff = 0.7,
                                    by = "p.adjust",
                                    select_fun = min)
write.table(as.data.frame(egoBP1), 'egoBP1.txt', sep = '\t', row.names = FALSE, quote = FALSE)
save(egoBP1,file='GOplot_input.Rdata')

##circle_dat()制作数据
load("../02_差异分析/DEG_limma_exprset51092.Rdata")
genelist <- nrDEG[tmp$SYMBOL,]
go <- egoBP1@result
colnames(go)
gocirc <- go[,c("ID","Description","p.adjust","geneID")]
colnames(gocirc)=c("ID","term","adj_pval","genes")
gocirc$genes <- gsub("/",",",gocirc$genes)
gocirc$category <- rep("BP",times = nrow(gocirc))
genelist$ID <- rownames(genelist)
##genelist <- merge(tmp,genelist,by.x='SYMBOL',by.y='ID')
##colnames(genelist)[colnames(genelist)=='ENTREZID'] <- 'ID'
circ <- circle_dat(gocirc, genelist)
head(circ)
##绘图
##原文：GOCircle(circ, nsub = c('GO:0009615', 'GO:0006955', 'GO:0006952','GO:0045087', 'GO:0034097', 'GO:0006954', 'GO:0043900', 'GO:0050670', 'GO:0070663', 'GO:0032944'))
##一些没有富集到，可能基因集和原文用的不一样
pdf("Fig3A_GOCircle.pdf",width = 10,height = 6)
GOCircle(circ, nsub = 10,label.size=4)
dev.off()

##展示基因与GO Terms关系的圈图 (GOChord())
head(circ)
##chord_dat ()将作图数据构建成GOChord() 要求的输入格式；一个二进制的关系矩阵，1表示基因属于该GO Term，0与之相反。
##选择感兴趣的基因
gochord <- gocirc[c(1:11),"ID"]
circ2 <- circ[circ$ID%in%gochord,]
genelist2 <- circ2[,c("genes","logFC")]
colnames(genelist2)[1]="ID"
genelist2 <- genelist2[!duplicated(genelist2$ID),]
chord <- chord_dat(data = circ2, genes = genelist2)
head(chord)
pdf("Fig3B_GOChord.pdf",width = 14,height = 14)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.2, gene.size = 3,border.size=0.01,process.label=10)
dev.off()
save(circ,circ2,chord,file='GOplot_output.Rdata')
