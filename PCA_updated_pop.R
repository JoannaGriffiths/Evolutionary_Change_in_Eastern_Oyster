


library("statmod")
library('DESeq2')
library('vegan')
library('ape')
library('raster')
library("edgeR")

install.packages("ggbiplot")
library("ggbiplot")
library("factoextra")

setwd("~/Documents/Oyster_exome/PCA")
setwd("~/LSU/Research/Oyster Exome Capture Experiment/PCA")

SL1_fam = read.delim("SL1-fam_exact_rc_pop_edit6", header = F)
colnames(SL1_fam) <- c('gene', 'SL3_top', 'SL3_bottom', 'SL4_top', 'SL4_bottom', 'SL1E_top', 'SL1E_bottom', 'SL1S_top', 'SL1S_bottom', 'SL3_freq', 'SL4_freq', 'SL1E_freq', 'SL1S_freq')
SL1_fam$SNP_ID <- row.names(SL1_fam)
SL1_fam$SL3_freq <- as.numeric(SL1_fam$SL3_freq)
SL1_fam$SL4_freq <- as.numeric(SL1_fam$SL4_freq)
SL1_fam$SL1E_freq <- as.numeric(SL1_fam$SL1E_freq)
SL1_fam$SL1S_freq <- as.numeric(SL1_fam$SL1S_freq)
SL1_fam$SNP_ID <- as.numeric(SL1_fam$SNP_ID)

SL1_fam.pca <- prcomp(SL1_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(SL1_fam.pca)
head(SL1_fam.pca$rotation)
SL1_fam_distance.matrix <- as.matrix(dist(SL1_fam.pca$rotation, method = "maximum", diag = T))

plot(SL1_fam.pca$rotation[,1], SL1_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")
#library('tdr')
#pca.dist <- pca2euclid(SL1_fam.pca)



SL2_fam = read.delim("SL2-fam_exact_rc_pop_edit6", header = F)
colnames(SL2_fam) <- c('gene', 'SL17_top', 'SL17_bottom', 'SL18_top', 'SL18_bottom', 'SL2E_top', 'SL2E_bottom', 'SL2S_top', 'SL2S_bottom', 'SL17_freq', 'SL18_freq', 'SL2E_freq', 'SL2S_freq')
SL2_fam$SNP_ID <- row.names(SL2_fam)
SL2_fam$SL17_freq <- as.numeric(SL2_fam$SL17_freq)
SL2_fam$SL18_freq <- as.numeric(SL2_fam$SL18_freq)
SL2_fam$SL2E_freq <- as.numeric(SL2_fam$SL2E_freq)
SL2_fam$SL2S_freq <- as.numeric(SL2_fam$SL2S_freq)
SL2_fam$SNP_ID <- as.numeric(SL2_fam$SNP_ID)

SL2_fam.pca <- prcomp(SL2_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(SL2_fam.pca)
head(SL2_fam.pca$rotation)

plot(SL2_fam.pca$rotation[,1], SL2_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")
library('tdr')
pca.dist <- pca2euclid(SL1_fam.pca)





SL3_fam = read.delim("SL3-fam_exact_rc_pop_edit7", header = F)
colnames(SL3_fam) <- c('gene', 'SL29_top', 'SL29_bottom', 'SL28_top', 'SL28_bottom', 'SL3E_top', 'SL3E_bottom', 'SL3S_top', 'SL3S_bottom', 'SL3NS_top', 'SL3NS_bottom', 'SL29_freq', 'SL28_freq', 'SL3E_freq', 'SL3S_freq', 'SL3NS_freq')
SL3_fam$SNP_ID <- row.names(SL3_fam)
SL3_fam$SL29_freq <- as.numeric(SL3_fam$SL29_freq)
SL3_fam$SL28_freq <- as.numeric(SL3_fam$SL28_freq)
SL3_fam$SL3E_freq <- as.numeric(SL3_fam$SL3E_freq)
SL3_fam$SL3S_freq <- as.numeric(SL3_fam$SL3S_freq)
SL3_fam$SL3NS_freq <- as.numeric(SL3_fam$SL3NS_freq)
SL3_fam$SNP_ID <- as.numeric(SL3_fam$SNP_ID)

SL3_fam.pca <- prcomp(SL3_fam[,c(12:16)], center = TRUE,scale. = TRUE)
summary(SL3_fam.pca)
head(SL3_fam.pca$rotation)

plot(SL3_fam.pca$rotation[,1], SL3_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red', 'purple'), xlab = "PCA 1", ylab = "PCA 2")



SL4_fam = read.delim("SL4-fam_exact_rc_pop_edit7", header = F)
colnames(SL4_fam) <- c('gene', 'SL56_top', 'SL56_bottom', 'SL58_top', 'SL58_bottom', 'SL4E_top', 'SL4E_bottom', 'SL4S_top', 'SL4S_bottom', 'SL4NS_top', 'SL4NS_bottom', 'SL56_freq', 'SL58_freq', 'SL4E_freq', 'SL4S_freq', 'SL4NS_freq')
SL4_fam$SNP_ID <- row.names(SL4_fam)
SL4_fam$SL56_freq <- as.numeric(SL4_fam$SL56_freq)
SL4_fam$SL58_freq <- as.numeric(SL4_fam$SL58_freq)
SL4_fam$SL4E_freq <- as.numeric(SL4_fam$SL4E_freq)
SL4_fam$SL4S_freq <- as.numeric(SL4_fam$SL4S_freq)
SL4_fam$SL4NS_freq <- as.numeric(SL4_fam$SL4NS_freq)
SL4_fam$SNP_ID <- as.numeric(SL4_fam$SNP_ID)

SL4_fam.pca <- prcomp(SL4_fam[,c(12:16)], center = TRUE,scale. = TRUE)
summary(SL4_fam.pca)
head(SL4_fam.pca$rotation)

plot(SL4_fam.pca$rotation[,1], SL4_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red', 'purple'), xlab = "PCA 1", ylab = "PCA 2")



AR2_fam = read.delim("AR2-fam_exact_rc_pop_edit6", header = F)
colnames(AR2_fam) <- c('gene', 'AR10_top', 'AR10_bottom', 'AR9_top', 'AR9_bottom', 'AR2E_top', 'AR2E_bottom', 'AR2S_top', 'AR2S_bottom', 'AR10_freq', 'AR9_freq', 'AR2E_freq', 'AR2S_freq')
AR2_fam$SNP_ID <- row.names(AR2_fam)
AR2_fam$AR10_freq <- as.numeric(AR2_fam$AR10_freq)
AR2_fam$AR9_freq <- as.numeric(AR2_fam$AR9_freq)
AR2_fam$AR2E_freq <- as.numeric(AR2_fam$AR2E_freq)
AR2_fam$AR2S_freq <- as.numeric(AR2_fam$AR2S_freq)
AR2_fam$SNP_ID <- as.numeric(AR2_fam$SNP_ID)

AR2_fam.pca <- prcomp(AR2_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(AR2_fam.pca)
head(AR2_fam.pca$rotation)

plot(AR2_fam.pca$rotation[,1], AR2_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



AR3_fam = read.delim("AR3-fam_exact_rc_pop_edit6", header = F)
colnames(AR3_fam) <- c('gene', 'AR55_top', 'AR55_bottom', 'AR13_top', 'AR13_bottom', 'AR3E_top', 'AR3E_bottom', 'AR3S_top', 'AR3S_bottom', 'AR55_freq', 'AR13_freq', 'AR3E_freq', 'AR3S_freq')
AR3_fam$SNP_ID <- row.names(AR3_fam)
AR3_fam$AR55_freq <- as.numeric(AR3_fam$AR55_freq)
AR3_fam$AR13_freq <- as.numeric(AR3_fam$AR13_freq)
AR3_fam$AR3E_freq <- as.numeric(AR3_fam$AR3E_freq)
AR3_fam$AR3S_freq <- as.numeric(AR3_fam$AR3S_freq)
AR3_fam$SNP_ID <- as.numeric(AR3_fam$SNP_ID)

AR3_fam.pca <- prcomp(AR3_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(AR3_fam.pca)
head(AR3_fam.pca$rotation)

plot(AR3_fam.pca$rotation[,1], AR3_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



AR4_fam = read.delim("AR4-fam_exact_rc_pop_edit6", header = F)
colnames(AR4_fam) <- c('gene', 'AR15_top', 'AR15_bottom', 'AR35_top', 'AR35_bottom', 'AR4E_top', 'AR4E_bottom', 'AR4S_top', 'AR4S_bottom', 'AR15_freq', 'AR35_freq', 'AR4E_freq', 'AR4S_freq')
AR4_fam$SNP_ID <- row.names(AR4_fam)
AR4_fam$AR15_freq <- as.numeric(AR4_fam$AR15_freq)
AR4_fam$AR35_freq <- as.numeric(AR4_fam$AR35_freq)
AR4_fam$AR4E_freq <- as.numeric(AR4_fam$AR4E_freq)
AR4_fam$AR4S_freq <- as.numeric(AR4_fam$AR4S_freq)
AR4_fam$SNP_ID <- as.numeric(AR4_fam$SNP_ID)

AR4_fam.pca <- prcomp(AR4_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(AR4_fam.pca)
head(AR4_fam.pca$rotation)

plot(AR4_fam.pca$rotation[,1], AR4_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



AR5_fam = read.delim("AR5-fam_exact_rc_pop_edit6", header = F)
colnames(AR5_fam) <- c('gene', 'AR1_top', 'AR1_bottom', 'AR29_top', 'AR29_bottom', 'AR5E_top', 'AR5E_bottom', 'AR5S_top', 'AR5S_bottom', 'AR1_freq', 'AR29_freq', 'AR5E_freq', 'AR5S_freq')
AR5_fam$SNP_ID <- row.names(AR5_fam)
AR5_fam$AR1_freq <- as.numeric(AR5_fam$AR1_freq)
AR5_fam$AR29_freq <- as.numeric(AR5_fam$AR29_freq)
AR5_fam$AR5E_freq <- as.numeric(AR5_fam$AR5E_freq)
AR5_fam$AR5S_freq <- as.numeric(AR5_fam$AR5S_freq)
AR5_fam$SNP_ID <- as.numeric(AR5_fam$SNP_ID)

AR5_fam.pca <- prcomp(AR5_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(AR5_fam.pca)
head(AR5_fam.pca$rotation)

plot(AR5_fam.pca$rotation[,1], AR5_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



VB2_fam = read.delim("VB2-fam_exact_rc_pop_edit6", header = F)
colnames(VB2_fam) <- c('gene', 'VB14_top', 'VB14_bottom', 'VB69_top', 'VB69_bottom', 'VB2E_top', 'VB2E_bottom', 'VB2S_top', 'VB2S_bottom', 'VB14_freq', 'VB69_freq', 'VB2E_freq', 'VB2S_freq')
VB2_fam$SNP_ID <- row.names(VB2_fam)
VB2_fam$VB14_freq <- as.numeric(VB2_fam$VB14_freq)
VB2_fam$VB69_freq <- as.numeric(VB2_fam$VB69_freq)
VB2_fam$VB2E_freq <- as.numeric(VB2_fam$VB2E_freq)
VB2_fam$VB2S_freq <- as.numeric(VB2_fam$VB2S_freq)
VB2_fam$SNP_ID <- as.numeric(VB2_fam$SNP_ID)

VB2_fam.pca <- prcomp(VB2_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(VB2_fam.pca)
head(VB2_fam.pca$rotation)

plot(VB2_fam.pca$rotation[,1], VB2_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



VB3R1_fam = read.delim("VB3_R1-fam_exact_rc_pop_edit6", header = F)
colnames(VB3R1_fam) <- c('gene', 'VB3_top', 'VB3_bottom', 'VB12_top', 'VB12_bottom', 'VB3E_top', 'VB3E_bottom', 'VB3S_top', 'VB3S_bottom', 'VB3_freq', 'VB12_freq', 'VB3E_freq', 'VB3S_freq')
VB3R1_fam$SNP_ID <- row.names(VB3R1_fam)
VB3R1_fam$VB3_freq <- as.numeric(VB3R1_fam$VB3_freq)
VB3R1_fam$VB12_freq <- as.numeric(VB3R1_fam$VB12_freq)
VB3R1_fam$VB3E_freq <- as.numeric(VB3R1_fam$VB3E_freq)
VB3R1_fam$VB3S_freq <- as.numeric(VB3R1_fam$VB3S_freq)
VB3R1_fam$SNP_ID <- as.numeric(VB3R1_fam$SNP_ID)

VB3R1_fam.pca <- prcomp(VB3R1_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(VB3R1_fam.pca)
head(VB3R1_fam.pca$rotation)

plot(VB3R1_fam.pca$rotation[,1], VB3R1_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")



VB3R2_fam = read.delim("VB3-R2fam_exact_rc_pop_edit6", header = F)
colnames(VB3R2_fam) <- c('gene', 'VB2_top', 'VB2_bottom', 'VB8_top', 'VB8_bottom', 'VB3E_top', 'VB3E_bottom', 'VB3S_top', 'VB3S_bottom', 'VB2_freq', 'VB8_freq', 'VB3E_freq', 'VB3S_freq')
VB3R2_fam$SNP_ID <- row.names(VB3R2_fam)
VB3R2_fam$VB2_freq <- as.numeric(VB3R2_fam$VB2_freq)
VB3R2_fam$VB8_freq <- as.numeric(VB3R2_fam$VB8_freq)
VB3R2_fam$VB3E_freq <- as.numeric(VB3R2_fam$VB3E_freq)
VB3R2_fam$VB3S_freq <- as.numeric(VB3R2_fam$VB3S_freq)
VB3R2_fam$SNP_ID <- as.numeric(VB3R2_fam$SNP_ID)

VB3R2_fam.pca <- prcomp(VB3R2_fam[,c(10:13)], center = TRUE,scale. = TRUE)
summary(VB3R2_fam.pca)
head(VB3R2_fam.pca$rotation)

plot(VB3R2_fam.pca$rotation[,1], VB3R2_fam.pca$rotation[,2], col=c('blue', 'orange', 'black', 'red'), xlab = "PCA 1", ylab = "PCA 2")




VB_pop = read.delim("VB-pop_exact_rc_pop_edit14", header = F)
VB_pop <- VB_pop[,c(26:37)]
colnames(VB_pop) <- c('VB14_freq', 'VB69_freq', 'VB2E_freq', 'VB2S_freq', 'VB3_freq', 'VB12_freq', 'VB3E_R1', 'VB3S_R1', 'VB2_freq', 'VB8_freq', 'VB3E_R2_freq', 'VB3S_R2_freq')
VB_pop$SNP_ID <- row.names(VB_pop)
VB_pop$SNP_ID <- as.numeric(VB_pop$SNP_ID)

VB_pop.pca <- prcomp(VB_pop[,c(1:12)], center = TRUE,scale. = TRUE)
summary(VB_pop.pca)
head(VB_pop.pca$rotation)
VB_pop_distance.matrix <- as.matrix(dist(VB_pop.pca$rotation, method = "maximum", diag = T))


quartz()
plot(VB_pop.pca$rotation[,1], VB_pop.pca$rotation[,2], 
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PCA 1", ylab = "PCA 2")

##########USE THIS
dds.pcoa=pcoa(vegdist(t((VB_pop[,1:12])),method="euclidean")/1000)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

quartz()
plot(scores[,1], scores[,2],  
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (31.87%)", ylab = "PC2 (26.83%)")

G1 <- pointDistance(dds.pcoa$vectors[,1:2],lonlat=F)
head(G1)

row.names(G1) <- row.names(dds.pcoa$vectors)
colnames(G1) <- row.names(dds.pcoa$vectors)
G1




################## 

AR_pop = read.delim("AR-pop_exact_rc_pop_edit18", header = F)
AR_pop <- AR_pop[,c(34:49)]
colnames(AR_pop) <- c('AR10_freq', 'AR9_freq', 'AR2E_freq', 'AR2S_freq', 'AR55_freq', 'AR13_freq', 'AR3E_freq', 'AR3S_freq', 'AR15_freq', 'AR35_freq', 'AR4E_freq', 'AR4S_freq', 'AR1_freq', 'AR29_freq', 'AR5E_freq', 'AR5S_freq')
AR_pop$SNP_ID <- row.names(AR_pop)
AR_pop$SNP_ID <- as.numeric(AR_pop$SNP_ID)

AR_pop.pca <- prcomp(AR_pop[,c(1:16)], center = TRUE,scale. = TRUE)
summary(AR_pop.pca)
head(AR_pop.pca$rotation)

AR_pop_distance.matrix <- as.matrix(dist(AR_pop.pca$rotation, method = "euclidean", diag = T))
#AR_pop_distance.matrix <- pointDistance(AR_pop.pca$rotation[,1:2],lonlat=TRUE)
rownames(AR_pop_distance.matrix) <- c('AR10_freq', 'AR9_freq', 'AR2E_freq', 'AR2S_freq', 'AR55_freq', 'AR13_freq', 'AR3E_freq', 'AR3S_freq', 'AR15_freq', 'AR35_freq', 'AR4E_freq', 'AR4S_freq', 'AR1_freq', 'AR29_freq', 'AR5E_freq', 'AR5S_freq')
colnames(AR_pop_distance.matrix) <- c('AR10_freq', 'AR9_freq', 'AR2E_freq', 'AR2S_freq', 'AR55_freq', 'AR13_freq', 'AR3E_freq', 'AR3S_freq', 'AR15_freq', 'AR35_freq', 'AR4E_freq', 'AR4S_freq', 'AR1_freq', 'AR29_freq', 'AR5E_freq', 'AR5S_freq')


##########USE THIS
dds.pcoa=pcoa(vegdist(t((AR_pop[,1:16])),method="euclidean")/1000)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

quartz()
plot(scores[,1], scores[,2],  
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'orange', 'orange', 'orange', 'orange', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (18.46%)", ylab = "PC2 (16.40%)")

G1 <- pointDistance(dds.pcoa$vectors[,1:2],lonlat=F)
head(G1)

row.names(G1) <- row.names(dds.pcoa$vectors)
colnames(G1) <- row.names(dds.pcoa$vectors)
G1
##################



quartz()
plot(AR_pop.pca$rotation[,1], AR_pop.pca$rotation[,2], 
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'orange', 'orange', 'orange', 'orange', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PCA 1", ylab = "PCA 2")
text(AR_pop ~colnames)


SL_pop = read.delim("SL-pop_exact_rc_pop_edit20", header = F)
SL_pop <- SL_pop[,c(38:55)]
colnames(SL_pop) <- c('SL3_freq', 'SL4_freq', 'SL1E_freq', 'SL1S_freq', 'SL17_freq', 'SL18_freq', 'SL2E_freq', 'SL2S_freq', 'SL29_freq', 'SL28_freq', 'SL3E_freq', 'SL3S_freq', 'SL3NS_freq', 'SL56_freq', 'SL58_freq', 'SL4E_freq', 'SL4S_freq', 'SL4NS_freq')
SL_pop <- SL_pop[c(1:13)]
SL_pop$SNP_ID <- row.names(SL_pop)
SL_pop$SNP_ID <- as.numeric(SL_pop$SNP_ID)

SL_pop.pca <- prcomp(SL_pop[,c(1:13)], center = TRUE,scale. = TRUE)
summary(SL_pop.pca)
head(SL_pop.pca$rotation)
SL_pop_distance.matrix <- as.matrix(dist(SL_pop.pca$rotation, method = "euclidean", diag = T))

windows()
plot(SL_pop.pca$rotation[,1], SL_pop.pca$rotation[,2], 
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red', 'red', 'orange', 'orange', 'orange', 'orange', 'orange'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 13),
     xlab = "PCA 1", ylab = "PCA 2")

##########USE THIS
dds.pcoa=pcoa(vegdist(t((SL_pop[,1:13])),method="euclidean")/1000)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

windows()
plot(scores[,1], scores[,2],  
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red', 'red', 'orange', 'orange', 'orange', 'orange', 'orange'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 13),
     xlab = "PC1 (28.12%)", ylab = "PC2 (21.88%)")

G1 <- pointDistance(dds.pcoa$vectors[,1:2],lonlat=F)
head(G1)

row.names(G1) <- row.names(dds.pcoa$vectors)
colnames(G1) <- row.names(dds.pcoa$vectors)
G1





################## NEW STUFF
setwd("~/Documents/Oyster_exome/PCA")
setwd("~/LSU/Research/Oyster Exome Capture Experiment/PCA")

all_pop = read.delim("all_exact_rc_edit", header = F) #round 1
colnames(all_pop) <- c('gene', 'SNP', 'AR10', 'AR9', 'AR2E', 'AR2S', 'AR55', 'AR13', 'AR3E', 'AR3S', 'AR15', 'AR35', 'AR4E', 'AR4S', 'AR1', 'AR29', 'AR5E', 'AR5S', 'SL3', 'SL4', 'SL1E', 'SL1S', 'SL17', 'SL18', 'SL2E', 'SL2S', 'SL29', 'SL28', 'SL3E', 'SL3S', 'SL3NS', 'SL56', 'SL58', 'SL4E', 'SL4S', 'SL4NS', 'VB14', 'VB69', 'VB2E', 'VB2S', 'VB3', 'VB12', 'VB3E_R1', 'VB3S_R1', 'VB2', 'VB8', 'VB3E_R2', 'VB3S_R2')
all_pop <- sapply(all_pop[,], as.character)
all_pop[,3:48] <- sapply(all_pop[,3:48], function(x) eval(parse(text=x))) #turns fractions into decimals

all_pop2 <- all_pop[,3:48]
all_pop2 <- all_pop2[,c(1:29,35:46)] #without missing data (SL4 family removed)
#all_pop2[,1:46] <- (sapply(all_pop2[,1:46], as.numeric)) #not entirely sure why this stopped working, but I think it has to do with character variables such as the decimal that prevents it from turning into numeric: https://stackoverflow.com/questions/2288485/how-to-convert-a-data-frame-column-to-numeric-type
sapply(all_pop2, class)

#run below code for round2 only
setwd("~/LSU/Research/Oyster Exome Capture Experiment/PCA")
#all_pop = read.delim("all_pop_exact_rc_edit2", header = T) #round 2
all_pop = read.delim("all_pop_exact_cov20_200_new_rc_edit2", header = T) #round 2
colnames(all_pop) <- c('gene', 'SNP', 'AR10', 'AR9', 'AR2E', 'AR2S', 'AR55', 'AR13', 'AR3E', 'AR3S', 'AR15', 'AR35', 'AR4E', 'AR4S', 'AR1', 'AR29', 'AR5E', 'AR5S', 'SL3', 'SL4', 'SL1E', 'SL1S', 'SL17', 'SL18', 'SL2E', 'SL2S', 'SL29', 'SL28', 'SL3E', 'SL3S', 'SL3NS', 'VB14', 'VB69', 'VB2E', 'VB2S', 'VB3', 'VB12', 'VB3E_R1', 'VB3S_R1', 'VB2', 'VB8', 'VB3E_R2', 'VB3S_R2')
all_pop <- sapply(all_pop[,], as.character)
all_pop[,3:43] <- sapply(all_pop[,3:43], function(x) eval(parse(text=x))) #turns fractions into decimals
all_pop2 <- all_pop[,3:43]
all_pop2[all_pop2 == "NaN"] <- 0

#round 1 and round 2 continue here
#cant figure out how to make code below trasnform whole data matrix without specifying each column
all_pop3<- transform(all_pop2, AR10 = as.numeric(AR10), AR9 = as.numeric(AR9), AR2E = as.numeric(AR2E), AR2S = as.numeric(AR2S), AR55 = as.numeric(AR55), AR13 = as.numeric(AR13), AR3E = as.numeric(AR3E), AR3S = as.numeric(AR3S), AR15 = as.numeric(AR15), AR35 = as.numeric(AR35), AR4E = as.numeric(AR4E), AR4S = as.numeric(AR4S), AR1 = as.numeric(AR1), AR29 = as.numeric(AR29), AR5E = as.numeric(AR5E), AR5S = as.numeric(AR5S), SL3 = as.numeric(SL3),
                     SL4 = as.numeric(SL4), SL1E = as.numeric(SL1E), SL1S = as.numeric(SL1S), SL17 = as.numeric(SL17), SL18 = as.numeric(SL18), SL2E = as.numeric(SL2E), SL2S = as.numeric(SL2S), SL29 = as.numeric(SL29), SL28 = as.numeric(SL28), SL3E = as.numeric(SL3E), SL3S = as.numeric(SL3S), SL3NS = as.numeric(SL3NS), 
                     VB14 = as.numeric(VB14), VB69 = as.numeric(VB69), VB2E = as.numeric(VB2E), VB2S = as.numeric(VB2S), VB3 = as.numeric(VB3), VB12 = as.numeric(VB12), VB3E_R1 = as.numeric(VB3E_R1), VB3S_R1 = as.numeric(VB3S_R1), VB2 = as.numeric(VB2), VB8 = as.numeric(VB8), VB3E_R2 = as.numeric(VB3E_R2), VB3S_R2 = as.numeric(VB3S_R2)
                        )
sapply(all_pop3, class) #checking its numeric

all_pop.pca <- prcomp(all_pop3[,c(1:41)], center = TRUE,scale. = TRUE)
all_pop.pca <- prcomp(all_pop3[c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,29,32,33,36,37,40,41)], center = TRUE,scale. = TRUE) #larvae only


#larvae only
larvae <- all_pop3[c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,29,32,33,36,37,40,41)] #larvae only
larvae.pca <- prcomp(larvae[,c(1:21)], center = TRUE,scale. = TRUE) #larvae only


#below section for trying to calculate distances which doesn't really work
summary(all_pop.pca) #get percent variance here
head(all_pop.pca$rotation) 
all_pop_distance.matrix <- as.matrix(dist(all_pop.pca$rotation, method = "euclidean", diag = T)) #not working since it gives me all the same distances

#Louisiana vs Texas populations
library("ellipse") #https://www.benjaminbell.co.uk/2018/02/principal-components-analysis-pca-in-r.html
# Get individuals (observations) as a matrix
tab <- matrix(c(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2]), ncol=2)
# Calculate correlations
LA_pop <- cor(tab[1:16,])
TX_pop <- cor(tab[17:41,])

larvae <- cor(tab[c(3:4,7:8,11:12,15:16,19:20,23:24,27:29,32:33,36:37,40:41),])
adults <- cor(tab[c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,30,31,34,35,38,39),])

windows()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2],  
     col=c('gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (16.1%)", ylab = "PC2 (8.9%)")
#TX vs LA pop, doesn't look great, don't include in results
polygon(ellipse(LA_pop/(max(abs(all_pop.pca$rotation))*100), centre=colMeans(tab[17:41,]), level=0.95), col=adjustcolor("skyblue", alpha.f=0.25), border="skyblue")
polygon(ellipse(TX_pop/(max(abs(all_pop.pca$rotation))*100), centre=colMeans(tab[1:16,]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold2")
##OR larvae vs adults, doesn't look great, don't include in results
polygon(ellipse(larvae/(max(abs(all_pop.pca$rotation))*100), centre=colMeans(tab[c(3:4,7:8,11:12,15:16,19:20,23:24,27:29,32:33,36:37,40:41),]), level=0.95), col=adjustcolor("skyblue", alpha.f=0.25), border="skyblue")
polygon(ellipse(adults/(max(abs(all_pop.pca$rotation))*100), centre=colMeans(tab[c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,30,31,34,35,38,39),]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold2")


#Louisiana vs Texas populations with true VB in blue
windows()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2],  
     col=c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (15.24%)", ylab = "PC2 (8.88%)")
text(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2], colnames( all_pop3 ), pos=3, cex=0.7)


#salinity exposure (also basically the same as year)
windows()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2],  
     col=c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (15.24%)", ylab = "PC2 (8.88%)")
text(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2], colnames( all_pop3 ), pos=3, cex=0.7)

#family
windows()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2],  
     col=c('black', 'black', 'black', 'black', 'yellow', 'yellow', 'yellow', 'yellow', 'red', 'red', 'red', 'red', 'green', 'green', 'green', 'green', 'orange', 'orange', 'orange', 'orange', 'aquamarine', 'aquamarine', 'aquamarine', 'aquamarine', 'purple', 'purple', 'purple', 'purple', 'purple', 'bisque1', 'bisque1', 'bisque1', 'bisque1', 'azure', 'azure', 'azure', 'azure', 'aliceblue', 'aliceblue', 'aliceblue', 'aliceblue'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19),
     xlab = "PC1 (16.1%)", ylab = "PC2 (8.9%)")
text(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2], colnames( all_pop3 ), pos=3, cex=0.7)
segments(x0 = all_pop.pca$rotation[3,1], y0 = all_pop.pca$rotation[3,2], x1 = all_pop.pca$rotation[4,1], y1 = all_pop.pca$rotation[4,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[7,1], y0 = all_pop.pca$rotation[7,2], x1 = all_pop.pca$rotation[8,1], y1 = all_pop.pca$rotation[8,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[11,1], y0 = all_pop.pca$rotation[11,2], x1 = all_pop.pca$rotation[12,1], y1 = all_pop.pca$rotation[12,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[15,1], y0 = all_pop.pca$rotation[15,2], x1 = all_pop.pca$rotation[16,1], y1 = all_pop.pca$rotation[16,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[19,1], y0 = all_pop.pca$rotation[19,2], x1 = all_pop.pca$rotation[20,1], y1 = all_pop.pca$rotation[20,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[23,1], y0 = all_pop.pca$rotation[23,2], x1 = all_pop.pca$rotation[24,1], y1 = all_pop.pca$rotation[24,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[27,1], y0 = all_pop.pca$rotation[27,2], x1 = all_pop.pca$rotation[28,1], y1 = all_pop.pca$rotation[28,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[27,1], y0 = all_pop.pca$rotation[27,2], x1 = all_pop.pca$rotation[29,1], y1 = all_pop.pca$rotation[29,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[32,1], y0 = all_pop.pca$rotation[32,2], x1 = all_pop.pca$rotation[33,1], y1 = all_pop.pca$rotation[33,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[36,1], y0 = all_pop.pca$rotation[36,2], x1 = all_pop.pca$rotation[37,1], y1 = all_pop.pca$rotation[37,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[40,1], y0 = all_pop.pca$rotation[40,2], x1 = all_pop.pca$rotation[41,1], y1 = all_pop.pca$rotation[41,2], col = "black") #SL1SE 



#larvae only, color is population
windows()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2],  
     col=c('gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2'), 
     pch = c(1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 13, 1, 19, 1, 19, 1, 19),
     xlab = "PC1 (16.1%)", ylab = "PC2 (8.9%)")
text(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2], colnames( larvae ), pos=3, cex=0.7)
segments(x0 = all_pop.pca$rotation[1,1], y0 = all_pop.pca$rotation[1,2], x1 = all_pop.pca$rotation[2,1], y1 = all_pop.pca$rotation[2,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[3,1], y0 = all_pop.pca$rotation[3,2], x1 = all_pop.pca$rotation[4,1], y1 = all_pop.pca$rotation[4,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[5,1], y0 = all_pop.pca$rotation[5,2], x1 = all_pop.pca$rotation[6,1], y1 = all_pop.pca$rotation[6,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[7,1], y0 = all_pop.pca$rotation[7,2], x1 = all_pop.pca$rotation[8,1], y1 = all_pop.pca$rotation[8,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[9,1], y0 = all_pop.pca$rotation[9,2], x1 = all_pop.pca$rotation[10,1], y1 = all_pop.pca$rotation[10,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[11,1], y0 = all_pop.pca$rotation[11,2], x1 = all_pop.pca$rotation[12,1], y1 = all_pop.pca$rotation[12,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[13,1], y0 = all_pop.pca$rotation[13,2], x1 = all_pop.pca$rotation[14,1], y1 = all_pop.pca$rotation[14,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[13,1], y0 = all_pop.pca$rotation[13,2], x1 = all_pop.pca$rotation[15,1], y1 = all_pop.pca$rotation[15,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[16,1], y0 = all_pop.pca$rotation[16,2], x1 = all_pop.pca$rotation[17,1], y1 = all_pop.pca$rotation[17,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[18,1], y0 = all_pop.pca$rotation[18,2], x1 = all_pop.pca$rotation[19,1], y1 = all_pop.pca$rotation[19,2], col = "black") #SL1SE 
segments(x0 = all_pop.pca$rotation[20,1], y0 = all_pop.pca$rotation[20,2], x1 = all_pop.pca$rotation[21,1], y1 = all_pop.pca$rotation[21,2], col = "black") #SL1SE 



##########USE THIS
library('vegan')
dds.pcoa=pcoa(vegdist(t((all_pop3)),method="manhattan")/1000) 
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

treatment_all <- read.delim("treatments_all.txt", header = T)
treatment_all$Salinity <- as.factor(treatment_all$Salinity)
treatment_all$Year <- as.factor(treatment_all$Year)
treatment_all$Population <- as.factor(treatment_all$Population)
treatment_all$Sig <- as.factor(treatment_all$Sig)
treatment_all$Time <- as.factor(treatment_all$Time)
treatment_all$LifeStage <- as.factor(treatment_all$LifeStage)
sapply(treatment_all, class)

grDevices::windows()
plot(scores[,1], scores[,2],  
     col=c('gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2'), 
     pch = c(0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1, 13, 0, 0, 19, 1, 0, 0, 19, 1, 0, 0, 19, 1),
     xlab = "PC1 (18.3%)", ylab = "PC2 (11.1%)", cex=1.5)
text(scores[,1], scores[,2], colnames( all_pop3 ), pos=3, cex=0.7)
segments(x0 = scores[3,1], y0 = scores[3,2], x1 = scores[4,1], y1 = scores[4,2], col = "black") #SL1SE 
segments(x0 = scores[7,1], y0 = scores[7,2], x1 = scores[8,1], y1 = scores[8,2], col = "black") #SL1SE 
segments(x0 = scores[11,1], y0 = scores[11,2], x1 = scores[12,1], y1 = scores[12,2], col = "black") #SL1SE 
segments(x0 = scores[15,1], y0 = scores[15,2], x1 = scores[16,1], y1 = scores[16,2], col = "black") #SL1SE 
segments(x0 = scores[19,1], y0 = scores[19,2], x1 = scores[20,1], y1 = scores[20,2], col = "black") #SL1SE 
segments(x0 = scores[23,1], y0 = scores[23,2], x1 = scores[24,1], y1 = scores[24,2], col = "black") #SL1SE 
segments(x0 = scores[27,1], y0 = scores[27,2], x1 = scores[28,1], y1 = scores[28,2], col = "black") #SL1SE 
segments(x0 = scores[27,1], y0 = scores[27,2], x1 = scores[29,1], y1 = scores[29,2], col = "black") #SL1SE 
segments(x0 = scores[32,1], y0 = scores[32,2], x1 = scores[33,1], y1 = scores[33,2], col = "black") #SL1SE 
segments(x0 = scores[36,1], y0 = scores[36,2], x1 = scores[37,1], y1 = scores[37,2], col = "black") #SL1SE 
segments(x0 = scores[40,1], y0 = scores[40,2], x1 = scores[41,1], y1 = scores[41,2], col = "black") #SL1SE 


adonis2(t(all_pop3)~Population, data=treatment_all, permutations = 1000000, method = "manhattan" ) # sig 0.000294 ***
adonis2(t(all_pop3)~LifeStage, data=treatment_all, permutations = 1000000, method = "manhattan" ) # sig 0.00769 **
adonis2(t(all_pop3)~Sig, data=treatment_all, permutations = 1000000, method = "manhattan" ) # sig 0.0158 *
adonis2(t(all_pop3)~Year, data=treatment_all, permutations = 1000000, method = "manhattan" )#  sig 1e-06 ***
adonis2(t(all_pop3)~Population+LifeStage+Sig+Year, data=treatment_all, permutations = 1000000, method = "manhattan" )#  sig 

#testing significance between group variances (need to do this because adonis will give sig p value even if groups overlap, but its sig b/c the variance is different)
#https://github.com/vegandevs/vegan/issues/233
#still getting this to work
short_all_pop3 <- head(all_pop3, n= 1000L)
pop_factor=as.character(treatment_all$Population)

input <- vegdist(t(all_pop3), method="manhattan", diag=FALSE, upper=TRUE, na.rm = TRUE)
mod <- betadisper(input, group = pop_factor, type = "median")
anova(mod)
#mod






##larvae only plot, color is population
treatment <- read.delim("treatments_larvae.txt", header = T)
treatment$Salinity <- as.factor(treatment$Salinity)
treatment$Year <- as.factor(treatment$Year)
treatment$Population <- as.factor(treatment$Population)
treatment$Sig <- as.factor(treatment$Sig)
treatment$Time <- as.factor(treatment$Time)
sapply(treatment, class)

larvae.dds.pcoa=pcoa(vegdist(t((larvae)),method="manhattan")/1000) 
larvae.scores=larvae.dds.pcoa$vectors
percent.larvae <- dds.pcoa$values$Eigenvalues
percent.larvae / sum(percent.larvae) #percent for each axes

grDevices::windows()
plot(larvae.scores[,1], larvae.scores[,2],  
     col=c('gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'gold', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2', 'skyblue2'), 
     pch = c(19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 13, 19, 1, 19, 1, 19, 1),
     xlab = "PC1 (18.3%)", ylab = "PC2 (11.1%)", cex=1.5)
text(larvae.scores[,1], larvae.scores[,2], colnames( larvae ), pos=3, cex=0.7)
segments(x0 = larvae.scores[1,1], y0 = larvae.scores[1,2], x1 = larvae.scores[2,1], y1 = larvae.scores[2,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[3,1], y0 = larvae.scores[3,2], x1 = larvae.scores[4,1], y1 = larvae.scores[4,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[5,1], y0 = larvae.scores[5,2], x1 = larvae.scores[6,1], y1 = larvae.scores[6,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[7,1], y0 = larvae.scores[7,2], x1 = larvae.scores[8,1], y1 = larvae.scores[8,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[9,1], y0 = larvae.scores[9,2], x1 = larvae.scores[10,1], y1 = larvae.scores[10,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[11,1], y0 = larvae.scores[11,2], x1 = larvae.scores[12,1], y1 = larvae.scores[12,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[13,1], y0 = larvae.scores[13,2], x1 = larvae.scores[14,1], y1 = larvae.scores[14,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[14,1], y0 = larvae.scores[14,2], x1 = larvae.scores[15,1], y1 = larvae.scores[15,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[16,1], y0 = larvae.scores[16,2], x1 = larvae.scores[17,1], y1 = larvae.scores[17,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[18,1], y0 = larvae.scores[18,2], x1 = larvae.scores[19,1], y1 = larvae.scores[19,2], col = "black") #SL1SE 
segments(x0 = larvae.scores[20,1], y0 = larvae.scores[20,2], x1 = larvae.scores[21,1], y1 = larvae.scores[21,2], col = "black") #SL1SE 


adonis2(t(larvae)~Population, data=treatment, permutations = 1000000, method = "manhattan" ) #significant, 0.004144 **
adonis2(t(larvae)~Time, data=treatment, permutations = 1000000, method = "manhattan" ) #not sig, 0.9943
adonis2(t(larvae)~Salinity, data=treatment, permutations = 1000000, method = "manhattan" ) # very sig, 2.5e-05 ***, but I think what matters is whether end point is sig or not
adonis2(t(larvae)~Salinity*Time, data=treatment, permutations = 1000000, method = "manhattan" ) # interaction not sig, 0.999
adonis2(t(larvae)~Sig, data=treatment, permutations = 1000000, method = "manhattan" ) # 0.08319 .
adonis2(t(larvae)~Sig*Time, data=treatment, permutations = 1000000, method = "manhattan" ) # not sig, 0.99
adonis2(t(larvae)~Year, data=treatment, permutations = 1000000, method = "manhattan" )# sig, 6e-05 ***
adonis2(t(larvae)~Population+Salinity*Time+Sig*Time+Year, data=treatment, permutations = 1000000, method = "manhattan" )# sig, 


pop_factor_larvae=as.character(treatment$Population)

input <- vegdist(t(larvae), method="manhattan", diag=FALSE, upper=TRUE, na.rm = TRUE)
mod <- betadisper(input, group = pop_factor_larvae, type = "median")
anova(mod)


G1 <- pointDistance(larvae.dds.pcoa$vectors[,1:2],lonlat=TRUE)
head(G1)

row.names(G1) <- row.names(larvae.dds.pcoa$vectors)
colnames(G1) <- row.names(larvae.dds.pcoa$vectors)
G1

distances <- read.delim("PCoA_larval_distances.txt", header=TRUE)
t.test(distances$Distance~distances$Sig) #not sig
t.test(distances$Distance~distances$Pop) #not sig
model1 = lm(Distance ~ Salinity, data = distances)
summary(aov(model1)) #not significant

##################

quartz()
plot(all_pop.pca$rotation[,1], all_pop.pca$rotation[,2], 
     col=c('blue', 'blue', 'blue', 'blue', 'black', 'black', 'black', 'black', 'red', 'red', 'red', 'red', 'red', 'orange', 'orange', 'orange', 'orange', 'orange'), 
     pch = c(0, 0, 1, 19, 0, 0, 1, 19, 0, 0, 1, 19, 13, 0, 0, 1, 19, 13),
     xlab = "PCA 1", ylab = "PCA 2")





######################
library(devtools)
install_github("vqv/ggbiplot")


vsd = assay(SL1_fam)
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes

#plot PC axis 1 and 2 for all data
plot(scores[,1], scores[,2], col=c("blue", "orange")[as.numeric(as.factor(factor2))], pch=c(0, 15, 1, 16)[as.numeric(as.factor(factor1))], xlab = "PC1 (42.47%)", ylab = "PC2 (17.22%)")
ordispider(scores,factor3, col=c("blue", "blue", "blue", "blue", "orange", "orange", "orange", "orange"))
