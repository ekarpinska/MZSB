matrix_1 <- read.csv(file = 'C:/Users/Ela/Downloads/mzsb/matrix.csv')
population <- read.csv(file ='C:/Users/Ela/Downloads/mzsb/igsr_populations.tsv', sep='\t')
(population$Superpopulation.name)
print(population)
print(matrix_1)
class(matrix_1)

#pca
rownames(matrix_1) <- matrix_1[,1]
dpca <-matrix_1[,2:5190]
pca <- prcomp(dpca)

dataset <- cbind(pca$x[,1], pca$x[,2], matrix_1[,5191])
data_ <- as.data.frame(dataset)
plot(data_[,1],data_[,2], col = as.factor(data_[,3]), xlab="PC1",ylab="PC2")
legend(10.5, 8.0, legend=c(unique(matrix_1[,5191])), 
       fill = as.factor(unique(matrix_1[,5191]))
)
pca
dataset <- cbind(pca$x[,1], pca$x[,2], matrix_1[,5191])
data_ <- as.data.frame(dataset)
(matrix_1[,5191])
population$Superpopulation.name

data_1 <- merge(dpca,(pca$x), by=0)
ggplot(data_1, aes(PC1, PC2, color=matrix_1[,5191])) + geom_point()

population$Superpopulation.code[population$Superpopulation.code == ""] <- "NI"

population$Superpopulation.code
plot(data_[,1],data_[,2], col = as.factor(population$Superpopulation.code), xlab="PC1",ylab="PC2")
legend("bottom", horiz=TRUE,xjust=1,yjust=0, 10.5, 8.0, legend=c(unique(population$Superpopulation.code)), 
       fill = as.factor(unique(population$Population.code))
)


#scree plot
library(tidyverse)
theme_set(theme_bw(24))


pca_variance <- pca$sdev^2
pca_calc <- round(pca_variance/sum(pca_variance)*100, 1)

barplot(pca_calc[1:5], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", names.arg=c("PC1", "PC2","PC3", "PC4", "PC5"))
pca_calc

#umap
iris.umap
install.packages("umap")

library(umap)
library("tidyverse")
library(umap)
BiocManager::install("M3C")
library(M3C)

iris.umap = umap(dpca)
#head(iris.umap$data[,1])
#plot(iris.umap, matrix_1[,5191])

plot((iris.umap$data[,1]),(iris.umap$data[,2]), col = as.factor(data_[,3]), xlab="x2",ylab="x1")

#boxplot
data2 <- data.frame(tsne_out$Y[,1], tsne_out$Y[,2], tsne_out$Y[,3], tsne_out$Y[,4], tsne_out$Y[,5])
barplot(data2, main="Scree Plot", xlab="t-sne Component", ylab="Percent Variation", names.arg=c("V1", "V2","V3", "V4", "V5"))

legend(angle=180, -30.0, 13.0, legend=c(unique(matrix_1[,5191])), 
        fill = as.factor(unique(matrix_1[,5191]))
)

plot((iris.umap$data[,1]),(iris.umap$data[,2]), col = as.factor(population$Superpopulation.code), xlab="x2",ylab="x1")


legend(angle=180, -2.0, 13.0, legend=c(unique(population$Superpopulation.code)), 
       fill = as.factor(unique(population$Population.code))
)
#umap_plot <- data.frame(x1 = (iris.umap$data[,1]), x2 = (iris.umap$data[,2]), col=population$Superpopulation.code)

#ggplot(umap_plot) + geom_point(aes(x=x1, y=x2, color=population$Superpopulation.code))

#population$Superpopulation.code

#as.matrix(iris.umap$data[,1]) -> iris.umap$data[,1]
#as.matrix(iris.umap$data[,2]) -> iris.umap$data[,2]
#nrow(iris.umap$data[,1])

#tsne
library(ggplot2)
library("tidyverse")
library(Rtsne)

tsne_out <- Rtsne(dpca) # TSNE

tsne_plot <- data.frame(tsn1 = tsne_out$Y[,1], tsn2 = tsne_out$Y[,2], col = matrix_1[,5191])
ggplot(tsne_plot) + geom_point(aes(x=tsn1, y=tsn2, color=col))

#boxplot
data2 <- data.frame(tsne_out$Y[,1], tsne_out$Y[,2], tsne_out$Y[,3], tsne_out$Y[,4], tsne_out$Y[,5])
barplot(data2, main="Scree Plot", xlab="t-sne Component", ylab="Percent Variation", names.arg=c("V1", "V2","V3", "V4", "V5"))


# sPCA
library("sparsepca")
sparse <- rspca(dpca, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
print(sparse)
summary(sparse)

sparse$scores


spca <- merge(dpca, (sparse$scores), by=0)
spca[,1:2]
library("tidyverse")
ggplot(spca, aes(V1, V2, color=matrix_1[,5191])) + geom_point()

#scree plot
x <- sparse$sdev^2
percent <- round(x/sum(x)*100, 1)
percent

barplot(per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", names.arg=c("V1", "V2","V3", "V4", "V5"))
sparse.var