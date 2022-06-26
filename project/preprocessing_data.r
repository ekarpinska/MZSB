library(vcfR)
install.packages("vcfR")
getwd()
vcf <- read.vcfR("ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")
vcf
my_data <- read.delim("phase1_integrated_calls.20101123.ALL.panel")
my_data
vcf
head(vcf)
getFIX(vcf)
vcf@gt
samples <-vcf@gt[,2:ncol(vcf@gt)]
ref <- getFIX(vcf)[,4]
ref
samples #counting allels
class(samples)
samples[samples =="0|0"] <- 0 #recessive homozygote
samples[samples =="0|1"] <- 1 
samples[samples =="1|0"] <- 1
samples[samples =="1|1"] <- 2
samples
ref_df <- as.data.frame(ref)
samples_df <- as.data.frame(samples)
ref_df[,1]
samples_df
row.names(samples_df)<-ref_df
samples_df = t(samples_df)
my_data[,1:2]
total <- merge(samples_df,my_data[,1:2],by=0)
write.csv(total,"C:/Users/Ela/Documents/matrix.csv", row.names = FALSE)

