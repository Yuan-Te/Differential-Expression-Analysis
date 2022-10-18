#Question 1
load('snps.RData')  # 53 samples x 1000 SNPs
summary(snps[,1:5])
nrow(snps)
ncol(snps)
head(snps)
colnames(snps)
rownames(snps)
class(snps)
snps

load('mat.gtex.RData')   #53 samples x 4 genes (expression)
summary(mat.gtex)
nrow(mat.gtex)
ncol(mat.gtex)
colnames(mat.gtex)
class(mat.gtex)

# Quality Control 
nSNPs <- ncol(snps) 
nSNPs
nSamples <- nrow(mat.gtex)
nSamples

samples.GenoData <- rownames(snps)
samples.ExprData <- rownames(mat.gtex)
setdiff(samples.GenoData, samples.ExprData)
setdiff(samples.ExprData, samples.GenoData)

length(intersect(samples.GenoData,samples.ExprData))
length(intersect(samples.GenoData,samples.ExprData)) == nSamples

hist(mat.gtex[,1], main="Gene expression profile: PIK3CA",
     xlab="Expression level")


df.mat.gtex <- data.frame(mat.gtex)
class(df.mat.gtex)
df.mat.gtex

library(ggplot2)
ggplot(df.mat.gtex, aes(PIK3CA))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
library(ggplot2)
ggplot(df.mat.gtex, aes(CDKN2A))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
library(ggplot2)
ggplot(df.mat.gtex, aes(TP53))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
library(ggplot2)
ggplot(df.mat.gtex, aes(SMAD4))+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

#Linear regression
df.mat.gtex <- data.frame(mat.gtex) #change class of mat.gtex from matrx array to dataframe
class(df.mat.gtex)
nrow(df.mat.gtex)
ncol(df.mat.gtex)

  #PIK3CA
pvalues.PIK3CA <- data.frame(MarkerName=colnames(snps),PIK3CA.pvalues=NA)

for (i in 1:ncol(snps)){ 
  mod <- lm(df.mat.gtex$PIK3CA~snps[,i])
  pvalues.PIK3CA[i,2] <- coef(summary(mod))[2,4]
}

head(pvalues.PIK3CA)
nrow(pvalues.PIK3CA)
ncol(pvalues.PIK3CA)

  #CDKN2A
pvalues.CDKN2A <- data.frame(MarkerName=colnames(snps),CDKN2A.pvalues=NA)

for (i in 1:ncol(snps)){ 
  mod <- lm(df.mat.gtex$CDKN2A~snps[,i])
  pvalues.CDKN2A[i,2] <- coef(summary(mod))[2,4]
}

head(pvalues.CDKN2A)
nrow(pvalues.CDKN2A)
ncol(pvalues.CDKN2A)

  #TP53
pvalues.TP53 <- data.frame(MarkerName=colnames(snps),TP53.pvalues=NA)

for (i in 1:ncol(snps)){ 
  mod <- lm(df.mat.gtex$TP53~snps[,i])
  pvalues.TP53[i,2] <- coef(summary(mod))[2,4]
}

head(pvalues.TP53)
nrow(pvalues.TP53)
ncol(pvalues.TP53)

  #SMAD4
pvalues.SMAD4 <- data.frame(MarkerName=colnames(snps),SMAD4.pvalues=NA)

for (i in 1:ncol(snps)){ 
  mod <- lm(df.mat.gtex$SMAD4~snps[,i])
  pvalues.SMAD4[i,2] <- coef(summary(mod))[2,4]
}

head(pvalues.SMAD4)
nrow(pvalues.SMAD4)
ncol(pvalues.SMAD4)

#Merge four linear regression
Total.pvalues <- data.frame(MarkerName=colnames(snps))
Total.pvalues

  #add PIK3CA
Total.pvalues <- merge(Total.pvalues, pvalues.PIK3CA, 
                                 by.x="MarkerName", by.y="MarkerName",
                                 all.x=FALSE, all.y=FALSE)

  #add CDKN2A  
Total.pvalues <- merge(Total.pvalues, pvalues.CDKN2A, 
                       by.x="MarkerName", by.y="MarkerName",
                       all.x=FALSE, all.y=FALSE)

  #add TP53  
Total.pvalues <- merge(Total.pvalues, pvalues.TP53, 
                       by.x="MarkerName", by.y="MarkerName",
                       all.x=FALSE, all.y=FALSE)

  #add SMAD4  
Total.pvalues <- merge(Total.pvalues, pvalues.SMAD4, 
                       by.x="MarkerName", by.y="MarkerName",
                       all.x=FALSE, all.y=FALSE)
ncol(Total.pvalues)
nrow(Total.pvalues)
head(Total.pvalues)


  #Using a cut-off of 10-8 to identify potential eQTLs for the 4 genes
sig.PIK3CA <- which(Total.pvalues$PIK3CA.pvalues<10^(-8))
length(sig.PIK3CA)
Total.pvalues[sig.PIK3CA,]

sig.CDKN2A <- which(Total.pvalues$CDKN2A.pvalues<10^(-8))
length(sig.CDKN2A)
Total.pvalues[sig.CDKN2A,]

sig.TP53 <- which(Total.pvalues$TP53.pvalues<10^(-8))
length(sig.TP53)
Total.pvalues[sig.TP53,]

sig.SMAD4 <- which(Total.pvalues$SMAD4.pvalues<10^(-8))
length(sig.SMAD4)
Total.pvalues[sig.SMAD4,]

# GGplot
 #Scatter plot
library(ggplot2)

scatter.PIK3CA.1 <- data.frame(PIK3CA=mat.gtex[,1], rs113505981=snps[,which(colnames(snps)=="rs113505981")])
scatter.PIK3CA.1
ggplot(scatter.PIK3CA.1, aes(rs113505981, PIK3CA))+
    geom_point()+
    geom_smooth(method=lm)


scatter.PIK3CA.2 <- data.frame(PIK3CA=mat.gtex[,1], rs4977264=snps[,which(colnames(snps)=="rs4977264")])
scatter.PIK3CA.2
ggplot(scatter.PIK3CA.2, aes(rs4977264, PIK3CA))+
  geom_point()+
  geom_smooth(method=lm)


scatter.SMAD4 <- data.frame(SMAD4=mat.gtex[,4], rs113505981=snps[,which(colnames(snps)=="rs113505981")])
scatter.SMAD4
ggplot(scatter.SMAD4, aes(rs113505981, SMAD4))+
  geom_point()+
  geom_smooth(method=lm)


 #Manhattan plot
plot(x = Teslovich.pvalues.SORT1$chrom_start,
     y=-log10(Teslovich.pvalues.SORT1$GC.Pvalue),
     pch='+',
     ylab='-log10(p) for biomarker',
     xlab='Position')





# Create a new covariate representing the tissue
list <- c('lung', 'pancreas', 'colon', 'oesophagus', 'brain', 'kidney')

Patient.PIK3CA <- data.frame(Patients=rownames(df.mat.gtex),
                             Expression=df.mat.gtex$PIK3CA,
                                    Tissue=NA)

setdiff(Patient.PIK3CA$Expression, df.mat.gtex$PIK3CA)
setdiff(df.mat.gtex$PIK3CA, Patient.PIK3CA$Expression)

  # Cases where PIK3CA expression is lower than 6 should only be found in pancreas, and cases with PIK3CA expression greater than 10 only in brain.
set.seed(4737346)
for (i in 1:nrow(Patient.PIK3CA)){
  if (Patient.PIK3CA[i,2]<6){
    Patient.PIK3CA[i,3] <- list[2]
  } else if (Patient.PIK3CA[i,2]>10){
    Patient.PIK3CA[i,3] <- list[5]
  } else {
    Patient.PIK3CA[i,3] <- list[sample(1:6, 1, replace=TRUE)]
  }
}

Patient.PIK3CA

  #Validation
a <- nrow(Patient.PIK3CA[Patient.PIK3CA$Expression>10,])
b <- nrow(Patient.PIK3CA[Patient.PIK3CA$Expression<6,])
c <- nrow(Patient.PIK3CA[Patient.PIK3CA$Expression<10 & Patient.PIK3CA$Expression>6,])
a+b+c
Patient.PIK3CA[Patient.PIK3CA$Expression>10,]
Patient.PIK3CA[Patient.PIK3CA$Expression<6,]
Patient.PIK3CA[Patient.PIK3CA$Expression<10 & Patient.PIK3CA$Expression>6,]




# Linear regression and a mixed effect model taking into account the tissue of origin

  #PIK3CA
pvalues.PIK3CA.tissue <- data.frame(MarkerName=colnames(snps),PIK3CA.pvalues=NA)

for (i in 1:ncol(snps)){ 
  mod <- lm(df.mat.gtex$PIK3CA~snps[,i]+Patient.PIK3CA[,3])
  pvalues.PIK3CA.tissue[i,2] <- coef(summary(mod))[2,4]
}

head(pvalues.PIK3CA.tissue)
nrow(pvalues.PIK3CA.tissue)
ncol(pvalues.PIK3CA.tissue)

sig.PIK3CA.tissue <- which(pvalues.PIK3CA.tissue$PIK3CA.pvalues<10^(-8))
length(sig.PIK3CA.tissue)
pvalues.PIK3CA.tissue[sig.PIK3CA.tissue,]    #from two set(without tissue) to 0 set(with tissue) 




#Question 2
load('tcga.coad.RData')

data('tcga.coad')

tcga.coad

q2 <- tcga.coad

summary(q2)
assay(q2)
nrow(q2)
colnames(q2) # 8 names of samples
rownames(q2) # 56602 genes
colSums(assay(q2))
head(assay(q2),3)
dim(q2)

str(metadata(rowRanges(q2))) #Use metadata to call database

attributes(q2) #Call out rowRanges & colData data

rowRanges(q2) #Call out data of rows(genes)

colData(q2) #Call out data of columns(samples)

  #Differential gene expression with DeSeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
q2$gender
q2$sample_type

head(assay(q2)) #等同head(assay(dds))，計算每組sample每個gene的reads
q2 <- q2[rowSums(assay(q2))>10, ]   #rowSums(counts(dds))將八組sample的reads加起來，只留下sum大於10的gene, 56602 rows只剩下26993 rows

head(assay(q2))
nrow(assay(q2))
ncol(assay(q2))  

colSums(assay(q2))


dds <- DESeqDataSet(q2, design = ~ 0+sample_type)  #BUG!!!!
colData(dds)

head(assay(dds))
attributes(dds)
colData(dds) 
rowRanges(dds)
class(dds) 

#normalization





