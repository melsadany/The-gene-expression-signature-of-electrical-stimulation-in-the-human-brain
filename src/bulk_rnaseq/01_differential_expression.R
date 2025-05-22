#Load librairies.
library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)

#Setup groups.
x <- as.factor(rep(c("baseline", "stimulated"), each=4))
names(x) <- colnames(Raw.Counts.Matrix)
groups <- matrix(data=c(1:4, 5:8), nrow=2, byrow=TRUE)

#Filter out undetected / lowly expressed genes.
filter <- apply(Raw.Counts.Matrix, 1, function(x) length(x[which(x>10)])>5)
filtered <- as.matrix(Raw.Counts.Matrix)[filter,]

#EDASeq normalization.
uq <- betweenLaneNormalization(filtered, which="upper")
#Observation of variation between samples
plotRLE(uq, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq, cex=1, cex.axis=1, cex.lab=1)

#EdgeR with RUV normalization.
controls <- rownames(uq)
s <- RUVs(uq, controls, k=2, groups)
plotRLE(s$normalizedCounts, outline=FALSE, las=3, ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(s$normalizedCounts, cex=1, cex.axis=1, cex.lab=1)

#designRUV
design <- model.matrix(~x + s$W)
design

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
topRsFC <- topTags(lrt, n=Inf)$table
