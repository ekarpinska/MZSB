library(devtools)
library(Biobase)
library(sva)
library(broom)
library(tidyverse)
library(data.table)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bladderbatch")
warnings()
#homework_1
## Attaching package: 'data.table'
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
## The following object is masked from 'package:purrr':
## 
##     transpose
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
save(bottomly.eset, file="bottomly.Rdata")

load(file="bottomly.Rdata")
ls()
## [1] "bottomly.eset" "con"           "setwd"
edata <- as.matrix(exprs(bottomly.eset))
dim(edata)
## [1] 36536    21
edata[1:5,1:5]
##                    SRX033480 SRX033488 SRX033481 SRX033489 SRX033482
## ENSMUSG00000000001       369       744       287       769       348
## ENSMUSG00000000003         0         0         0         0         0
## ENSMUSG00000000028         0         1         0         1         1
## ENSMUSG00000000031         0         0         0         0         0
## ENSMUSG00000000037         0         1         1         5         0
edata <- edata[rowMeans(edata) > 10, ]
edata <- log2(as.matrix(edata) + 1)
edata <- t(scale(t(edata), scale=FALSE, center=TRUE))
svd.out <- svd(edata)

PC = data.table(svd.out$v,pData(bottomly.eset))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(lane.number)))
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(experiment.number)))
print(cor(pData(bottomly.eset)$experiment.number, svd.out$v[,1], method="spearman"))

print(cor(pData(bottomly.eset)$experiment.number, svd.out$v[,2], method="spearman"))
print(cor(pData(bottomly.eset)$lane.number, svd.out$v[,1], method="spearman"))
print(cor(pData(bottomly.eset)$lane.number, svd.out$v[,2], method="spearman"))

library(bladderbatch)
data(bladderdata)

# sample info
pheno = pData(bladderEset)
# expression data
edata = exprs(bladderEset)

dim(pheno)
dim(edata)
edata[1:5,1:10]
head(pheno) 
sumna <- apply(edata, 1, function(x) sum(is.na(x)))
row.variances <- apply(edata, 1, function(x) var(x))
row.means <- apply(edata, 1, function(x) mean(x))
plot(row.variances, row.means, pch=19, main="Mean vs. Variance relationship")
edata <- edata[row.variances < 6,]
edata.log <- log2(edata)
edata.scaled <- t(scale(t(edata.log), scale=TRUE, center=TRUE))
edata.centered <- t(scale(t(edata.log), scale=FALSE, center=TRUE))

svd.centered.out <- svd(edata.centered)
svd.centered.plot <- data.table(svd.centered.out$v[,1:10], pheno)

svd.scaled.out <- svd(edata.scaled)
svd.scaled.plot <- data.table(svd.scaled.out$v[,1:10], pheno)
ggplot(svd.centered.plot) + geom_point(aes(x=V1, y=V2, col=as.factor(batch)))
ggplot(svd.centered.plot) + geom_point(aes(x=V1, y=V2, col=as.factor(cancer)))
ggplot(svd.scaled.plot) + geom_point(aes(x=V1, y=V2, col=as.factor(batch)))
ggplot(svd.scaled.plot) + geom_point(aes(x=V1, y=V2, col=as.factor(cancer)))

pheno_batch <- pheno[order(pheno$batch),]

pheno_batch <- pheno_batch[, c("batch", "outcome")]

table(pheno_batch)

#homework_2
pheno$cancer = relevel(pheno$cancer, ref="Normal")
mod = lm(edata[1,] ~ as.factor(pheno$cancer) + as.factor(pheno$batch))
print(mod)
pheno$cancer = relevel(pheno$cancer, ref="Normal")
mod = lm(t(edata) ~ as.factor(pheno$cancer) + as.factor(pheno$batch))
names(mod)
dim(mod$coefficients)
rownames(mod$coefficients)

# library "broom" clean up the outputs of LM
# now, we can use ggplot2 to plot various aspects of LM
library(broom)
mod_tidy <- tidy(mod)
library(RColorBrewer)
ggplot(mod_tidy) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")
# however, the previous line of code make a histogram of all coefficients.
# what we need to do is to find estimates of particular regression terms.
mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")
ggplot(mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")
# how about the p-values?
ggplot(mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")
library(sva)
batch = pheno$batch
combat_edata = ComBat(dat=edata, batch=pheno$batch, mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)
class(combat_edata)
dim(combat_edata)
combat_edata[1:10,1:10]
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)

fit <- kmeans(edata, 5)
aggregate(edata,by=list(fit$cluster),FUN=mean)
data <- data.frame(edata, fit$cluster)
head(data)
kmeans <- data %>% arrange(fit.cluster)
head(kmeans)
a<-as.matrix(kmeans)
head(a)
png("bladder.png",height=700,width=700)
heatmap.2(a,
          main = "Bladder Cancer Data Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "col",
          Colv=FALSE)
dev.off()
png("bladder.png",height=700,width=700)
fit <- kmeans(combat_edata, 5)
aggregate(combat_edata,by=list(fit$cluster),FUN=mean)

data <- data.frame(combat_edata, fit$cluster)
kmeans <- data %>% arrange(fit.cluster)
a<-as.matrix(kmeans)
getwd
png("bladder_combat.png",height=700,width=700)
heatmap.2(a,
          main = "Bladder Cancer Data Cleaned by ComBat", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "col")
dev.off()
svd.out.combat <- svd(combat_edata)
svd.combat.plot <- data.table(svd.out.combat$v[,1:10], pheno)

ggplot(svd.combat.plot) + geom_point(aes(x=V1, y=V2, col=as.factor(batch)))

modcombat = lm(t(combat_edata) ~ as.factor(pheno$cancer))

# library "broom" clean up the outputs of LM
# now, we can use ggplot2 to plot various aspects of LM
library(broom)
modcombat_tidy <- tidy(modcombat)

# histogram of estimates of particular regression terms.
ggplot(modcombat_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")
a

#homework_3
head(pheno)
head(edata)

edata.Pearson = cor(edata, method = 'pearson')
head(edata.Pearson)
head(edata)
edata

study.ident = rownames(edata.Pearson)

# set labels
label = list()
for (i in 1:length(study.ident)){
  label[i] <- paste(study.ident[i], pheno[,2][i], pheno[,3][i])
}

rownames(edata.Pearson) <- label
colnames(edata.Pearson) <- label

my_palette <- colorRampPalette(c("blue", "white", 'yellow', "darkred"))(n=299)
png("pearson.png",height=700,width=700)

heatmap.2(edata.Pearson,
          main = "Correlations given labels of phenotype data", 
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          scale = 'none',
          Colv=FALSE)
dev.off()
getwd()


#homework_4
new_pheno <- pData(bottomly.eset)
new_pheno
linear_edata = lm(t(edata) ~ as.factor(new_pheno$strain) + as.factor(new_pheno$experiment.number))
linear_edata_tidy <- tidy(linear_edata)

library(sva)
library(broom)
batch = new_pheno$batch
combat_edata = ComBat(dat=edata, batch=new_pheno$experiment.number,
                      mod=model.matrix(~1, data=new_pheno), par.prior=TRUE, prior.plots=TRUE)
modcombat <- lm(t(combat_edata) ~ as.factor(new_pheno$strain))


# library "broom" clean up the outputs of LM
# now, we can use ggplot2 to plot various aspects of LM

modcombat_tidy <- tidy(modcombat)

# histogram of estimates of particular regression terms.
ggplot(modcombat_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J")) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")
new_pheno
# different way of looking at coefficients from many different models
# the vertical line indicates the zero coefficient.
ggplot(modcombat_tidy, aes(estimate, term)) +
  geom_point() +
  geom_vline(xintercept = 0)

# filter : choose ROWS 
# select : choose COLS
est_compare <- tibble(
  LinearModel = linear_edata_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("estimate") %>% unlist,
  ComBat = modcombat_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("estimate") %>% unlist)

ggplot(est_compare, aes(x=LinearModel, y=ComBat)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()



pvalues.compare <- tibble(
  Linear_Model_ = linear_edata_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("p.value") %>% unlist,
  ComBat_ = modcombat_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("p.value") %>% unlist)
ggplot(pvalues.compare, aes(x= Linear_Model_)) + geom_histogram()
ggplot(pvalues.compare, aes(x= ComBat_), add=T) + geom_histogram()


#homework_5

mod = model.matrix(~as.factor(strain), data=new_pheno)

# permutation procedure from Buja and Eyuboglu 1992

num.sv(edata,mod,method="be")

# asymptotic approach from Leek 2011 Biometrics.
num.sv(edata,mod,method="leek")

mod = model.matrix(~as.factor(strain),data=new_pheno)
mod0 = model.matrix(~1, data=new_pheno)
sva_output = sva(edata, mod, mod0, n.sv=num.sv(edata,mod,method="leek"))

head(sva_output$sv)

# summary shows how the batches are related to SV1 and SV2 separately.
summary(lm(sva_output$sv ~ new_pheno$experiment.number))

sva_batch <- tibble(SV1=sva_output$sv[,1],
                    SV2=sva_output$sv[,2],
                    SV3=sva_output$sv[,3],
                    SV4=sva_output$sv[,4],
                    batch=as.factor(pheno$batch),
                    cancer=as.factor(pheno$cancer),
                    outcome=as.factor(pheno$outcome))

ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=batch))
# Add the surrogate variables to the model matrix
modsva = lm(t(edata) ~ as.factor(new_pheno$strain) + sva_output$sv)
modsva_tidy <- tidy(modsva)

est_compare <- tibble(
  LinearModel = linear_edata_tidy %>% filter(term == "as.factor(new_pheno$cancer)DBA/2J") %>% select("estimate") %>% unlist,
  
  ComBat = modcombat_tidy %>% filter(term == "as.factor(new_pheno$cancer)DBA/2J") %>% select("estimate") %>% unlist,
  
  SVA = modsva_tidy %>% filter(term == "as.factor(new_pheno$cancer)DBA/2J") %>% select("estimate") %>% unlist)

ggplot(est_compare, aes(x=LinearModel, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

ggplot(est_compare, aes(x=ComBat, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

ggplot(modsva_tidy %>% filter(term == "as.factor(new_pheno$cancer)DBA/2J")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")

pvalues <- tibble(
  LinearModel = linear_edata_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("p.value") %>% unlist,
  ComBat = modcombat_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("p.value") %>% unlist,
  SVA = modsva_tidy %>% filter(term == "as.factor(new_pheno$strain)DBA/2J") %>% select("p.value") %>% unlist)

pvalues_gather <- gather(pvalues)
ggplot(pvalues_gather, aes(x=value)) + geom_histogram() + facet_wrap(~key)

