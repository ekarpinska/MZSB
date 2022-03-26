#problem1

library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)
library(tidyverse)
library(data.table)
library(corpcor)

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
print(con)
load(file=con)
close(con)
summary(bottomly.eset)
class(bottomly.eset)
save(bottomly.eset, file="bottomly.Rdata")
load(file="bottomly.Rdata")

edata <- as.matrix(exprs(bottomly.eset))
dim(edata)
edata <- log2(as.matrix(edata) + 1)
edata <- edata[rowMeans(edata) > 10, ]
library(RColorBrewer)
library(gplots)
my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 299)

png("bottomly_heatmap_raw.png",height=700,width=700)
heatmap.2(edata,
          main = "Bottomly et al. Raw", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",    # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)
dev.off()
png("bottomly_heatmap_clustered.png",height=700,width=700)

wss <- (nrow(edata)-1)*sum(apply(edata,2,var))

fit <- kmeans(edata, 5)

aggregate(edata,by=list(fit$cluster),FUN=mean)

data <- data.frame(edata, fit$cluster)
data_matrix<-as.matrix(data)

data <- data_matrix[,1:21]
heatmap.2((data),
          main = "Bottomly et al. Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="col",     # only draw a row dendrogram
          scale = "col"
          )
dev.off()

#problem2

edata <- t(scale(t(edata), scale=FALSE, center=TRUE))
svd.out <- svd(edata)
names(svd.out)

print(paste("Dimension of left singular vectors:", dim(svd.out$u)))

print(paste("Length of singular values:",length(svd.out$d)))

print(paste("Dimension of right singular vectors:",dim(svd.out$v)))

par(mfrow=c(1,2))
plot(svd.out$d, pch=20, ylab="Singular values")
plot(svd.out$d^2/sum(svd.out$d^2)*100, pch=20, ylab="% variance explained")


plot(1:ncol(edata),svd.out$v[,1],pch=20)
library(data.table)
PC<-data.table(svd.out$v,pData(bottomly.eset))
head(PC)

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(strain)))

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(lane.number)))

ggplot(PC)+geom_point(aes(x=V1,y=V2,col=as.factor(experiment.number)))

ggplot(PC)+geom_point(aes(x=V3,y=V5,col=as.factor(strain)))

# problem3

PC<-data.table(svd.out$u,pData(bottomly.eset)) #left singular
ggplot(PC) + geom_boxplot(aes(x=as.factor(strain), y=V1))
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75))
ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V2),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V2))
head(PC)
ggplot(PC, aes(loadings, amount)) + geom_jitter(aes(x="V1", y=V1, col=as.factor(strain)))+  geom_jitter(aes(x="V2", y=V2, col=as.factor(strain)))

#problem 4

ggplot(PC) + geom_violin(aes(x=as.factor(strain), y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x=as.factor(strain), y=V1))
head(PC[, c("V1", "V2","V3","V4","V5")])

table <- (PC[, c("V1", "V2","V3","V4","V5")])

new_data <- as.data.frame(table)
head(new_data)
ggplot(new_data, aes(Loadings, Amount))+ geom_violin(aes(x="V1", y=V1),draw_quantiles = c(0.25, 0.5, 0.75))+ geom_violin(aes(x="V1", y=V1),draw_quantiles = c(0.25, 0.5, 0.75)) + geom_violin(aes(x="V2", y=V2),draw_quantiles = c(0.25, 0.5, 0.75))+ geom_violin(aes(x="V3", y=V3),draw_quantiles = c(0.25, 0.5, 0.75))+ geom_violin(aes(x="V4", y=V4),draw_quantiles = c(0.25, 0.5, 0.75))+ geom_violin(aes(x="V5", y=V5),draw_quantiles = c(0.25, 0.5, 0.75))
