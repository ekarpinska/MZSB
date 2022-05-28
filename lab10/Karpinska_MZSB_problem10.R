install.packages("Seurat")
install.packages("Rtsne")
library(dplyr)
library(Seurat)
library(patchwork)
library(Rtsne)
library(tidyverse)

pbmc.data <- Read10X(data.dir = "C:/Users/Ela/OneDrive/Documents/filtered_matrices_mex/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = 'percent.mt')

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = 'umap')

#Homework 1
pbmc_em <- Embeddings(object = pbmc, reduction = "pca")[,1:15]

wss <- (nrow(pbmc_em)-1)*sum(apply(pbmc_em,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(pbmc_em,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(pbmc_em, 10)

t_sne <- Rtsne(pbmc_em,initial_dims=50, perplexity=30)
t_sne
t_snee = data.frame(t_sne$Y)
t_snee
t_snee$cluster <- k_means$cluster
ggplot(t_snee) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")
k_means$cluster

#homework2
#clu <- list(pbmc_em[k_means$cluster == 1,], pbmc_em[k_means$cluster == 2,], pbmc_em[k_means$cluster == 3,], pbmc_em[k_means$cluster == 4,], pbmc_em[k_means$cluster == 5,], pbmc_em[k_means$cluster == 6,], pbmc_em[k_means$cluster == 7,], pbmc_em[k_means$cluster == 8,], pbmc_em[k_means$cluster == 9,], pbmc_em[k_means$cluster == 10,])
clu1 <- pbmc_em[k_means$cluster == 1,]
clu2 <- pbmc_em[k_means$cluster == 2,]
clu3 <- pbmc_em[k_means$cluster == 3,]
clu4 <- pbmc_em[k_means$cluster == 4,]
clu5 <- pbmc_em[k_means$cluster == 5,]
clu6 <- pbmc_em[k_means$cluster == 6,]
clu7 <- pbmc_em[k_means$cluster == 7,]
clu8 <- pbmc_em[k_means$cluster == 8,]
clu9 <- pbmc_em[k_means$cluster == 9,]
clu10 <- pbmc_em[k_means$cluster == 10,]
#clu[1]
wss <- (nrow(clu1)-1)*sum(apply(clu1,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu1),
                                              centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu1, 8)
t_sne <- Rtsne(as.matrix(clu1), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[2]
wss <- (nrow(clu2)-1)*sum(apply(clu2,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu2),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu2, 2)
t_sne <- Rtsne(as.matrix(clu2), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[3]
wss <- (nrow(clu3)-1)*sum(apply(clu3,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu3),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu3, 6)
t_sne <- Rtsne(as.matrix(clu3), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[4]
wss <- (nrow(clu4)-1)*sum(apply(clu4,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu4),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu4, 5)
t_sne <- Rtsne(as.matrix(clu4), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[5]
wss <- (nrow(clu5)-1)*sum(apply(clu5,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu5),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu5, 4)
t_sne <- Rtsne(as.matrix(clu5), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[6]
wss <- (nrow(clu6)-1)*sum(apply(clu6,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu6),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu6, 6)
t_sne <- Rtsne(as.matrix(clu6), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[7]
wss <- (nrow(clu7)-1)*sum(apply(clu7,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu7),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu7, 5)
t_sne <- Rtsne(as.matrix(clu7), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[8]
wss <- (nrow(clu8)-1)*sum(apply(clu8,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu8),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu8, 8)
t_sne <- Rtsne(as.matrix(clu8), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[9]
wss <- (nrow(clu9)-1)*sum(apply(clu9,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu9),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu9, 7)
t_sne <- Rtsne(as.matrix(clu9), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

#clu[10]
wss <- (nrow(clu10)-1)*sum(apply(clu10,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(clu10),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(clu10, 4)
t_sne <- Rtsne(as.matrix(clu10), perplexity = 30)
t_sne = data.frame(t_sne$Y)

t_sne$cluster <- k_means$cluster
ggplot(t_sne) + geom_point(aes(x=X1, y=X2, col = as.factor(cluster)))  + labs(color = "cluster")

