


library("statmod")
library('DESeq2')
library('vegan')
library('ape')
library('raster')
library("edgeR")

install.packages("ggbiplot")
library("ggbiplot")
library("factoextra")


setwd("~/LSU/Research/Oyster Exome Capture Experiment/PCA")

all_pop = read.delim("all_pop_exact_cov20_200_new_rc_edit2", header = T) 
colnames(all_pop) <- c('gene', 'SNP', 'AR10', 'AR9', 'AR2E', 'AR2S', 'AR55', 'AR13', 'AR3E', 'AR3S', 'AR15', 'AR35', 'AR4E', 'AR4S', 'AR1', 'AR29', 'AR5E', 'AR5S', 'SL3', 'SL4', 'SL1E', 'SL1S', 'SL17', 'SL18', 'SL2E', 'SL2S', 'SL29', 'SL28', 'SL3E', 'SL3S', 'SL3NS', 'VB14', 'VB69', 'VB2E', 'VB2S', 'VB3', 'VB12', 'VB3E_R1', 'VB3S_R1', 'VB2', 'VB8', 'VB3E_R2', 'VB3S_R2')
all_pop <- sapply(all_pop[,], as.character)
all_pop[,3:43] <- sapply(all_pop[,3:43], function(x) eval(parse(text=x))) #turns fractions into decimals
all_pop2 <- all_pop[,3:43]
all_pop2[all_pop2 == "NaN"] <- 0

all_pop3<- transform(all_pop2, AR10 = as.numeric(AR10), AR9 = as.numeric(AR9), AR2E = as.numeric(AR2E), AR2S = as.numeric(AR2S), AR55 = as.numeric(AR55), AR13 = as.numeric(AR13), AR3E = as.numeric(AR3E), AR3S = as.numeric(AR3S), AR15 = as.numeric(AR15), AR35 = as.numeric(AR35), AR4E = as.numeric(AR4E), AR4S = as.numeric(AR4S), AR1 = as.numeric(AR1), AR29 = as.numeric(AR29), AR5E = as.numeric(AR5E), AR5S = as.numeric(AR5S), SL3 = as.numeric(SL3),
                     SL4 = as.numeric(SL4), SL1E = as.numeric(SL1E), SL1S = as.numeric(SL1S), SL17 = as.numeric(SL17), SL18 = as.numeric(SL18), SL2E = as.numeric(SL2E), SL2S = as.numeric(SL2S), SL29 = as.numeric(SL29), SL28 = as.numeric(SL28), SL3E = as.numeric(SL3E), SL3S = as.numeric(SL3S), SL3NS = as.numeric(SL3NS), 
                     VB14 = as.numeric(VB14), VB69 = as.numeric(VB69), VB2E = as.numeric(VB2E), VB2S = as.numeric(VB2S), VB3 = as.numeric(VB3), VB12 = as.numeric(VB12), VB3E_R1 = as.numeric(VB3E_R1), VB3S_R1 = as.numeric(VB3S_R1), VB2 = as.numeric(VB2), VB8 = as.numeric(VB8), VB3E_R2 = as.numeric(VB3E_R2), VB3S_R2 = as.numeric(VB3S_R2)
                        )
sapply(all_pop3, class) #checking its numeric


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

#Calculating distances in PCoA plot
G1 <- pointDistance(larvae.dds.pcoa$vectors[,1:2],lonlat=TRUE)
head(G1)

row.names(G1) <- row.names(larvae.dds.pcoa$vectors)
colnames(G1) <- row.names(larvae.dds.pcoa$vectors)
G1

#determing if distances are correlated with treatment factors
distances <- read.delim("PCoA_larval_distances.txt", header=TRUE)
distances$Sig <- as.numeric(distances$Sig)
summary(aov(lm(Distance ~ Sig, data = distances))) #not sig
summary(aov(lm(Distance ~ Pop, data = distances))) #not sig
summary(aov(lm(Distance ~ Salinity, data = distances))) #not significant

