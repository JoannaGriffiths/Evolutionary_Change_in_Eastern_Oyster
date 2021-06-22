
library("dplyr")
setwd("~/LSU/Research/Oyster Exome Capture Experiment")

###
annot_final2 <- read.delim("~/exome_desktop/Supplemental_TableS1", header = T)


###



#####AR2
AR2 <- read.delim("AR2-SE_sig_edit5.fet", header = F)
colnames(AR2) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

AR2$SNP_pos <- AR2$Start + AR2$SNP
AR2$Punlog <- (1/(10^(AR2$p_value)))
AR2$padj <- p.adjust(AR2$Punlog, method = 'fdr')

AR2_sig <- subset(AR2, padj < 0.05)
write.table(AR2_sig, file = "AR2_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

AR2_count <- AR2_sig %>% count(Gene) #count number of significant SNPs per gene

mean(AR2_count$n)
densities.qtiles <- AR2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(AR2_count$n)
AR2_count3 <- subset(AR2_count, n>=3)
colnames(AR2_count3) <- c("Gene", "AR2_n")

AR2_annot <- merge(AR2_count3, annot_final2, by="Gene")


#####AR3 (none sig)
AR3 <- read.delim("AR3-SE_sig_edit5.fet", header = F)
colnames(AR3) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

AR3$SNP_pos <- AR3$Start + AR3$SNP
AR3$Punlog <- (1/(10^(AR3$p_value)))
AR3$padj <- p.adjust(AR3$Punlog, method = 'fdr')

AR3_sig <- subset(AR3, padj < 0.05)
write.table(AR3_sig, file = "AR3_sig_SNPs_forSeqMonk", row.names=FALSE, quote=FALSE)

AR3_count <- AR3_sig %>% count(Gene) #count number of significant SNPs per gene

mean(AR3_count$n)
densities.qtiles <- AR3_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(AR3_count$n)
AR3_count3 <- subset(AR3_count, n>=3)
colnames(AR3_count3) <- c("Gene", "AR3_n")

AR3_annot <- merge(AR3_count3, annot_final2, by="Gene")


#####AR4
AR4 <- read.delim("AR4-SE_sig_edit5.fet", header = F)
colnames(AR4) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

AR4$SNP_pos <- AR4$Start + AR4$SNP
AR4$Punlog <- (1/(10^(AR4$p_value)))
AR4$padj <- p.adjust(AR4$Punlog, method = 'fdr')

AR4_sig <- subset(AR4, padj < 0.05)
write.table(AR4_sig, file = "AR4_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

AR4_count <- AR4_sig %>% count(Gene) #count number of significant SNPs per gene

mean(AR4_count$n)
densities.qtiles <- AR4_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(AR4_count$n)
AR4_count3 <- subset(AR4_count, n>=3)
colnames(AR4_count3) <- c("Gene", "AR4_n")

AR4_annot <- merge(AR4_count3, annot_final2, by="Gene")


#####AR5
AR5 <- read.delim("AR5-SE_sig_edit5.fet", header = F)
colnames(AR5) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

AR5$SNP_pos <- AR5$Start + AR5$SNP
AR5$Punlog <- (1/(10^(AR5$p_value)))
AR5$padj <- p.adjust(AR5$Punlog, method = 'fdr')

AR5_sig <- subset(AR5, padj < 0.05)
write.table(AR5_sig, file = "AR5_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

AR5_count <- AR5_sig %>% count(Gene) #count number of significant SNPs per gene

mean(AR5_count$n)
densities.qtiles <- AR5_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(AR5_count$n)
AR5_count3 <- subset(AR5_count, n>=3)
colnames(AR5_count3) <- c("Gene", "AR5_n")

AR5_annot <- merge(AR5_count3, annot_final2, by="Gene")


#####SL1
SL1 <- read.delim("SL1-SE_sig_edit5.fet", header = F)
colnames(SL1) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

SL1$SNP_pos <- SL1$Start + SL1$SNP
SL1$Punlog <- (1/(10^(SL1$p_value)))
SL1$padj <- p.adjust(SL1$Punlog, method = 'fdr')

SL1_sig <- subset(SL1, padj < 0.05)
write.table(SL1_sig, file = "SL1_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

SL1_count <- SL1_sig %>% count(Gene) #count number of significant SNPs per gene

mean(SL1_count$n)
densities.qtiles <- SL1_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(SL1_count$n)
SL1_count3 <- subset(SL1_count, n>=3)
colnames(SL1_count3) <- c("Gene", "SL1_n")

SL1_annot <- merge(SL1_count3, annot_final2, by="Gene")


#####SL2 (none significant)
SL2 <- read.delim("SL2-SE_sig_edit5.fet", header = F)
colnames(SL2) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

SL2$SNP_pos <- SL2$Start + SL2$SNP
SL2$Punlog <- (1/(10^(SL2$p_value)))
SL2$padj <- p.adjust(SL2$Punlog, method = 'fdr')

SL2_sig <- subset(SL2, padj < 0.05)
write.csv2(SL2_sig, file = "SL2_sig_SNPs_forSeqMonk", row.names=FALSE, quote=FALSE)

SL2_count <- SL2_sig %>% count(Gene) #count number of significant SNPs per gene

mean(SL2_count$n)
densities.qtiles <- SL2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(SL2_count$n)
SL2_count3 <- subset(SL2_count, n>=3)
colnames(SL2_count3) <- c("Gene", "SL2_n")

SL2_annot <- merge(SL2_count3, annot_final2, by="Gene")


#####SL3 (none significant)
SL3 <- read.delim("SL3-SE_sig_edit5.fet", header = F)
colnames(SL3) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

SL3$SNP_pos <- SL3$Start + SL3$SNP
SL3$Punlog <- (1/(10^(SL3$p_value)))
SL3$padj <- p.adjust(SL3$Punlog, method = 'fdr')

SL3_sig <- subset(SL3, padj < 0.05)


SL3_count <- SL3_sig %>% count(Gene) #count number of significant SNPs per gene

mean(SL3_count$n)
densities.qtiles <- SL3_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(SL3_count$n)
SL3_count3 <- subset(SL3_count, n>=3)
colnames(SL3_count3) <- c("Gene", "SL3_n")

SL3_annot <- merge(SL3_count3, annot_final2, by="Gene")


#####SLNS3 (none significant)
SLNS3 <- read.delim("SLNS3-SE_sig_edit5.fet", header = F)
colnames(SLNS3) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

SLNS3$SNP_pos <- SLNS3$Start + SLNS3$SNP
SLNS3$Punlog <- (1/(10^(SLNS3$p_value)))
SLNS3$padj <- p.adjust(SLNS3$Punlog, method = 'fdr')

SLNS3_sig <- subset(SLNS3, padj < 0.05)

SLNS3_count <- SLNS3_sig %>% count(Gene) #count number of significant SNPs per gene

mean(SLNS3_count$n)
densities.qtiles <- SLNS3_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(SLNS3_count$n)
SLNS3_count3 <- subset(SLNS3_count, n>=10)
colnames(SLNS3_count3) <- c("Gene", "SLNS3_n")

SLNS3_annot <- merge(SLNS3_count3, annot_final2, by="Gene")


#####VB2
VB2 <- read.delim("VB2-SE_sig_edit5.fet", header = F)
colnames(VB2) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

VB2$SNP_pos <- VB2$Start + VB2$SNP
VB2$Punlog <- (1/(10^(VB2$p_value)))
VB2$padj <- p.adjust(VB2$Punlog, method = 'fdr')

VB2_sig <- subset(VB2, padj < 0.05)
write.table(VB2_sig, file = "VB2_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

VB2_count <- VB2_sig %>% count(Gene) #count number of significant SNPs per gene

mean(VB2_count$n)
densities.qtiles <- VB2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(VB2_count$n)
VB2_count3 <- subset(VB2_count, n>=3)
colnames(VB2_count3) <- c("Gene", "VB2_n")

VB2_annot <- merge(VB2_count3, annot_final2, by="Gene")


#####VB3
VB3 <- read.delim("VB3-SE_sig_edit5.fet", header = F)
colnames(VB3) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

VB3$SNP_pos <- VB3$Start + VB3$SNP
VB3$Punlog <- (1/(10^(VB3$p_value)))
VB3$padj <- p.adjust(VB3$Punlog, method = 'fdr')

VB3_sig <- subset(VB3, padj < 0.05)
write.table(VB3_sig, file = "VB3_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

VB3_count <- VB3_sig %>% count(Gene) #count number of significant SNPs per gene

mean(VB3_count$n)
densities.qtiles <- VB3_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(VB3_count$n)
VB3_count3 <- subset(VB3_count, n>=3)
colnames(VB3_count3) <- c("Gene", "VB3_n")

VB3_annot <- merge(VB3_count3, annot_final2, by="Gene")


#####VB3_R2
VB3_R2 <- read.delim("VB3_R2-SE_sig_edit5.fet", header = F)
colnames(VB3_R2) <- cbind("Gene","Species","CHR","Start","End","SNP","IDK1","IDK2", "coverage", "p_value") #adding column names to dataframe

VB3_R2$SNP_pos <- VB3_R2$Start + VB3_R2$SNP
VB3_R2$Punlog <- (1/(10^(VB3_R2$p_value)))
VB3_R2$padj <- p.adjust(VB3_R2$Punlog, method = 'fdr')

VB3_R2_sig <- subset(VB3_R2, padj < 0.05)
write.table(VB3_R2_sig, file = "VB3_R2_sig_SNPs", row.names=FALSE, quote=FALSE, sep = "\t")

VB3_R2_count <- VB3_R2_sig %>% count(Gene) #count number of significant SNPs per gene

mean(VB3_R2_count$n)
densities.qtiles <- VB3_R2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

windows()
hist(VB3_R2_count$n)
VB3_R2_count3 <- subset(VB3_R2_count, n>=3)
colnames(VB3_R2_count3) <- c("Gene", "VB3_R2_n")

VB3_R2_annot <- merge(VB3_R2_count3, annot_final2, by="Gene")


##############################chisq test for source enrichment
AR2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(0,1,1, 1,0))
chisq.test(AR2_df) #not sig 0.73


AR4_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(0,2,1, 1,1))
chisq.test(AR4_df) #not sig 0.89


AR5_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(3,3,4, 1,4))
chisq.test(AR5_df) #not sig 0.75


SL1_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(0,0,2,0,0))
chisq.test(SL1_df) #not sig 0.07


VB2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(1,0,0, 0,2))
chisq.test(VB2_df) #not sig 0.31

VB3_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(1,1,0,0,0))
chisq.test(VB3_df) #not sig 0.52


VB3R2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(5,5,8,5,8))
chisq.test(VB3R2_df) #not sig 0.33


all_combined_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(10,12,16, 8,15))
chisq.test(all_combined_df) #not sig 0.75



AR2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,1,0,0,1,0,1,0,1,0,0,0))
chisq.test(AR2_df) #not sig 0.61


AR4_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,1,1,0,0,1,0,0,0,0,1,0))
chisq.test(AR4_df) #not sig 0.82


AR5_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(2,0,3,0,0,1,0,1,2,0,1,1))
chisq.test(AR5_df) #not sig 0.74


SL1_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(1,0,0,0,0,0,0,0,1,0,0,0))
chisq.test(SL1_df) #not sig 0.72


VB2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,0,0,0,1,0,0,0,0,0,1,1))
chisq.test(VB2_df) #not sig 0.91

VB3R1_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,0,0,0,0,0,0,0,0,0,0,1))
chisq.test(VB3R1_df) #not sig 0.78

VB3R2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(5,4,6,0,3,2,1,0,2,0,3,1))
chisq.test(VB3R2_df) #not sig 0.93





########## Hypergeometric test
'''
x, q: vector of quantiles representing the number of white balls drawn without replacement 
          from an urn which contains both black and white balls.
m: the number of white balls in the urn.
n:the number of black balls in the urn.
k:the number of balls drawn from the urn.

q: number of shared significant genes among all crosses
m: total number of significant genes for a particular cross? between two crosses?
n: total number of possible genes minus total number of significant genes
k: total number significant genes (added up among all crosses?)

https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper 
https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
'''

tot_sig_genes <- merge(AR2_count3, AR4_count3, by="Gene", all = T)
tot_sig_genes2 <- merge(tot_sig_genes, AR5_count3, by="Gene", all = T)
tot_sig_genes3 <- merge(tot_sig_genes2, VB2_count3, by="Gene", all = T)
tot_sig_genes4 <- merge(tot_sig_genes3, SL1_count3, by="Gene", all = T)
tot_sig_genes5 <- merge(tot_sig_genes4, VB3_R2_count3, by="Gene", all = T)
tot_sig_genes6 <- merge(tot_sig_genes5, VB3_count3, by="Gene", all = T)

#43 unique genes total


#AR2: 4sig
#AR3: 0sig
#AR4: 4sig
#AR5: 13sig
#SL1: 2
#SL2: 0
#VB2: 3
#VB3R1: 1
#VB3R2: 30

#phyper((number of significant genes in this cross that overlap with other crosses),
#(number of genes significant in other crosses), 
#(number of genes not significant in other crosses), 
#(number of genes significant in this cross that don't overlap with any other cross)) 

#AR2 (TX1) 3 overlap, 43-(4-3), 152-42, 4-3
1-(phyper(3, 42, 110, 1)) #not sig 0

#AR4 (TX3) 2 overlap, 43-(4-2), 152-41, 4-2
1-(phyper(2, 41, 111, 2)) #sig at 0

#AR5 (TX4) 7 overlap, 43-(13-7), 152-37, 13-7
1-(phyper(7, 37, 115, 6)) #sig at 0

#VB2 (LA1) 1 overlap, 43-(3-1), 152-41, 3-1
1-(phyper(1, 41, 111, 2)) #not sig at 0.07

#VB3R2 (LA3) 11 overlap, 43-(30-11), 152-24, 30-11
1-(phyper(11, 24, 128, 19)) #not sig at 3.720706e-07
