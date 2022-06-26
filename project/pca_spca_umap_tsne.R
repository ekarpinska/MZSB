library(dplyr)
matrix_1 <- read.csv(file = 'C:/Users/Ela/Documents/Ela_mzsb/matrix.csv')
population <- read.csv(file ='C:/Users/Ela/Documents/Ela_mzsb/igsr_populations.tsv', sep='\t')

df <- data.frame(population$Population.code, population$Superpopulation.name)
colnames(df) <- c("Population.code","Superpopulation.name")
data_merge1 <- merge(matrix_1, df, by = "Population.code") 
data_merge1[,5192]
Countries <- matrix_1[,5191]
########## PCA countries
matrix_1[,5191]
rownames(matrix_1) <- matrix_1[,1]
matrix_1[1:6,1:6]
dpca <-matrix_1[,2:5190]
pca <- prcomp(dpca)

dataset <- cbind(pca$x[,1], pca$x[,2], data_merge1[,5192])
data_ <- as.data.frame(dataset)
pca$x[1:2,]
dataset

pca
dataset <- cbind(pca$x[,1], pca$x[,2], Countries)
data_ <- as.data.frame(dataset)
(matrix_1[,5191])
population$Superpopulation.name
library("tidyverse")
data_1 <- merge(dpca,(pca$x), by=0)
ggplot(data_1, aes(PC1, PC2, color=Countries)) + geom_point()

###### PCA super
Superpolulations <- data_merge1[,5192]
dataset <- cbind(pca$x[,1], pca$x[,2], Superpolulations)
data_ <- as.data.frame(dataset)

data_1 <- merge(dpca,(pca$x), by=0)
ggplot(data_1, aes(PC1, PC2, color=Superpolulations)) + geom_point()


## Scree plot

theme_set(theme_bw(24))

pca$sdev
pca_variance <- pca$sdev^2
pca_calc <- round(pca_variance/sum(pca_variance)*100, 1)
pca_calc
barplot(pca_calc[1:5], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", names.arg=c("PC1", "PC2","PC3", "PC4", "PC5"))

#umap countries

install.packages("umap")

library(umap)
dft<-umap.defaults
dft$n_neighbors = 15
dft

a<-umap(dpca,config=dft, method = c("naive", "umap-learn"),
        preserve.seed = TRUE)

umap_plot <- data.frame(x1 = (a$layout[,1]), x2 = (a$layout[,2]), Countries=Countries)

ggplot(umap_plot) + geom_point(aes(x=x1, y=x2, color=Countries))

##umap super

umap_plot <- data.frame(x1 = (a$layout[,1]), x2 = (a$layout[,2]), Superpolulations=Superpolulations)

ggplot(umap_plot) + geom_point(aes(x=x1, y=x2, color=Superpolulations))


#tSNE countries
library(Rtsne)

tsne_out <- Rtsne(dpca)
tsne_plot <- data.frame(tsn1 = tsne_out$Y[,1], tsn2 = tsne_out$Y[,2], Countries = Countries)
ggplot(tsne_plot) + geom_point(aes(x=tsn1, y=tsn2, color=Countries))

#tSNE Superpopulations

tsne_plot <- data.frame(tsn1 = tsne_out$Y[,1], tsn2 = tsne_out$Y[,2], Superpolulations = Superpolulations)
ggplot(tsne_plot) + geom_point(aes(x=tsn1, y=tsn2, color=Superpolulations))


# sPCA Countries
library("sparsepca")
sparse <- rspca(dpca, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
spca <- merge(dpca, (sparse$scores), by=0)
spca
ggplot(spca, aes(V1, V2, color=Countries)) + geom_point()

# sPCA Superpopulations
ggplot(spca, aes(V1, V2, color=Superpolulations)) + geom_point()
sparse$objective
##sPCA Scree plot
var_expl<-summary(sparse)[3,]*100
var_expl
barplot(var_expl, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", names.arg=c("PC1", "PC2","PC3", "PC4", "PC5"))


###### Transposition
#transpose matrix
matrix_2 <- t(matrix_1)
print (matrix_2)
colnames(matrix_2) <- matrix_2[1,]
matrix_2[1:6,1:6]
as.data.frame(matrix_2[2:6,1:6])
transpose <- matrix_2[2:5190,1:1092]
t1=data.frame(transpose)
t1[5189,1092]
t1 <- as.numeric(t1)
i <- c(1, 5189)  
t1[i,] <- apply(t1[i,], 1,            # Specify own function within apply
                function(x) as.numeric(as.character(x)))
sapply(t1[1:6,1:6], class)
t2 <- data.frame(sapply(t1, function(x) as.numeric(as.character(x))))
pca_snp <- prcomp(t2)
pca_snp$x[1:5,1:5]
colnames(t1)
rownames(t1)<- t1[,1]
matrix_2[5191,1:6]
dataset <- cbind(pca_snp$x[,1], pca_snp$x[,2], data_merge1[,5192])
data_ <- as.data.frame(dataset)
pca_snp$x[1:2,]
dataset
plot(data_[,1],data_[,2], col = as.factor(data_[,3]), xlab="PC1",ylab="PC2")
legend(50, 22.0, legend=c(unique(data_merge1[,5192])), 
       fill = as.factor(unique(data_merge1[,5192]))
)

###################
population
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


####################spca transpose
sparse <- rspca(t2, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
print(sparse)
summary(sparse)

sparse$scores
dataset <- cbind(spca[,1094], spca[,1095], matrix_2[5191,])
data_ <- as.data.frame(dataset)
pca$x[1:2,]
dataset
plot(data_[,1],data_[,2], col = as.factor(data_[,3]), xlab="PC1",ylab="PC2")
legend(40, 3.0, legend=c(unique(matrix_1[,5191])), 
       fill = as.factor(unique(matrix_1[,5191]))
)


#finish of transpose analysis
                        
#kmeans
wss <- (nrow(dpca)-1)*sum(apply(dpca,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(dpca),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(dpca, 4)

umap_plot <- data.frame(x1 = (iris.umap$data[,1]), x2 = (iris.umap$data[,2]), col=k_means$cluster)

ggplot(umap_plot) + geom_point(aes(x=x1, y=x2, color=k_means$cluster))
                        
                        
#tsne k means


wss <- (nrow(dpca)-1)*sum(apply(dpca,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(dpca),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(dpca, 3)

#umap kmeans
wss <- (nrow(dpca)-1)*sum(apply(dpca,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(as.matrix(dpca),
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
k_means <- kmeans(dpca, 4)

umap_plot <- data.frame(x1 = (a.umap$data[,1]), x2 = (a.umap$data[,2]), col=k_means$cluster)

ggplot(umap_plot) + geom_point(aes(x=x1, y=x2, color=k_means$cluster))
