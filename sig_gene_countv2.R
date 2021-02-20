
library("dplyr")

setwd("C:/Users/joann/OneDrive/Documents/exome_desktop/fet_sliding_window")
setwd("C:/Users/joann/Documents/oyster_exome/")
annot <- read.delim("C:/Users/joann/Desktop/oyster_exome/uniq_gene_list2", header = T)

setwd("~/oyster_exome")
annot <- read.delim("uniq_gene_list2", header = T)
setwd("~/exome_desktop/fet_sliding_window")

AR2 <- read.delim("AR2-newpvalue_Manhattan")
AR2_sig <- subset(AR2, padj < 0.05)

AR2_count <- AR2_sig %>% count(Gene) #count number of significant SNPs per gene

mean(AR2_count$n)
densities.qtiles <- AR2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(AR2_count$n)
AR2_count3 <- subset(AR2_count, n>=3)
colnames(AR2_count3) <- c("Gene", "AR2_n")

AR2_annot <- merge(AR2_count3, annot, by="Gene")



AR3 <- read.delim("AR3-newpvalue_Manhattan")
AR3_sig <- subset(AR3, padj < 0.05)

AR3_count <- AR3_sig %>% count(Gene)

mean(AR3_count$n)
densities.qtiles <- AR3_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(AR3_count$n)
AR3_count3 <- subset(AR3_count, n>=3) #most is 2

AR3_annot <- merge(AR3_count3, annot, by="Gene")



AR4 <- read.delim("AR4-newpvalue_Manhattan")
AR4_sig <- subset(AR4, padj < 0.05)

AR4_count <- AR4_sig %>% count(Gene)

mean(AR4_count$n)
densities.qtiles <- AR4_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(AR4_count$n)
AR4_count3 <- subset(AR4_count, n>=3)
colnames(AR4_count3) <- c("Gene", "AR4_n")

AR4_annot <- merge(AR4_count3, annot, by="Gene")




AR5 <- read.delim("AR5-newpvalue_Manhattan")
AR5_sig <- subset(AR5, padj < 0.05)

AR5_count <- AR5_sig %>% count(Gene)

mean(AR5_count$n)
densities.qtiles <- AR5_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(AR5_count$n)
AR5_count3 <- subset(AR5_count, n>=3)
colnames(AR5_count3) <- c("Gene", "AR5_n")

AR5_annot <- merge(AR5_count3, annot, by="Gene")




SL1 <- read.delim("SL1-newpvalue_Manhattan")
SL1_sig <- subset(SL1, padj < 0.05)

SL1_count <- SL1_sig %>% count(Gene)

mean(SL1_count$n)
densities.qtiles <- SL1_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(SL1_count$n)
SL1_count3 <- subset(SL1_count, n>=3)
colnames(SL1_count3) <- c("Gene", "SL1_n")

SL1_annot <- merge(SL1_count3, annot, by="Gene")



SL2 <- read.delim("SL2-newpvalue_Manhattan")
SL2_sig <- subset(SL2, padj < 0.05)

SL2_count <- SL2_sig %>% count(Gene)

mean(SL2_count$n)
densities.qtiles <- SL2_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(SL2_count$n)
SL2_count3 <- subset(SL2_count, n>=3)
colnames(SL2_count3) <- c("Gene", "SL2_n")

SL2_annot <- merge(SL2_count3, annot, by="Gene")



VB3R1 <- read.delim("VB3R1-newpvalue_Manhattan")
VB3R1_sig <- subset(VB3R1, padj < 0.05)

VB3R1_count <- VB3R1_sig %>% count(Gene)

mean(VB3R1_count$n)
densities.qtiles <- SL1_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(VB3R1_count$n) #nothing past count 2


VB3R2 <- read.delim("VB3R2-newpvalue_Manhattan")
VB3R2_sig <- subset(VB3R2, padj < 0.05)

VB3R2_count <- VB3R2_sig %>% count(gene)
colnames(VB3R2_count) <- c("Gene", "VB3R2_n")

mean(VB3R2_count$n)
densities.qtiles <- SL1_count %>%
  summarise(q05 = quantile(n, 0.025, na.rm=TRUE),
            q50 = quantile(n, 0.5, na.rm=TRUE),
            q95 = quantile(n, 0.975, na.rm=TRUE))

hist(VB3R2_count$n)
VB3R2_count3 <- subset(VB3R2_count, n>=3)
colnames(VB3R2_count3) <- c("Gene", "VB3R2_n")

VB3R2_annot <- merge(VB3R2_count3, annot, by="Gene")

##############################chisq test
AR2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(3,2,3, 3,1))
chisq.test(AR2_df) #not sig 0.55


AR4_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(2,1,1, 1,0))
chisq.test(AR4_df) #not sig 0.50


AR5_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(2,2,2, 0,3))
chisq.test(AR5_df) #not sig 0.67


SL1_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(0,1,1,0,0))
chisq.test(SL1_df) #not sig 0.64


SL2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(1,1,0, 1,0))
chisq.test(SL2_df) #not sig 0.65


VB3R2_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(10,14,13,4,18))
chisq.test(VB3R2_df) #not sig 0.33


all_combined_df <- data.frame(
  row.names=c("Fst", "exp", "she", "meng", "filler"),
  all=c(26,47,33, 30,42),
  cross=c(18,21,20, 9,11))
chisq.test(all_combined_df) #not sig 0.14



AR2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(1,2,2,1,0,1,1,0,1,1,0,1))
chisq.test(AR2_df) #not sig 0.30


AR4_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,1,1,0,0,1,0,0,0,0,0,1))
chisq.test(AR4_df) #not sig 0.84


AR5_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,3,1,0,0,1,1,1,0,0,0,1))
chisq.test(AR5_df) #not sig 0.38


SL1_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(1,0,0,0,0,0,0,0,1,0,0,0))
chisq.test(SL1_df) #not sig 0.72


SL2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(0,0,1,0,0,0,0,0,0,0,0,1))
chisq.test(SL2_df) #not sig 0.93


VB3R2_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(5,7,11,0,4,3,3,1,4,0,4,7))
chisq.test(VB3R2_df) #not sig 0.93


all_combined_df <- data.frame(
  row.names=c("1chemical_def", "2immune", "3ion", "4ELF", "5FAA_meta", "6FA_hydro", "7ubiq", "8proteolysis", "9ROS", "10RNA_poly", "11FAA_protein", "12uncharacterized"),
  all=c(18,22,27,3,22,7,5,4,13,1,15,18),
  cross=c(7,13,16,1,4,5,5,2,6,1,4,11))
chisq.test(all_combined_df) #not sig 0.69



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
tot_sig_genes3 <- merge(tot_sig_genes2, SL2_count3, by="Gene", all = T)
tot_sig_genes4 <- merge(tot_sig_genes3, SL1_count3, by="Gene", all = T)
tot_sig_genes5 <- merge(tot_sig_genes4, VB3R2_count3, by="Gene", all = T)



#AR2: 11sig
#AR3: 0sig
#AR4: 4sig
#AR5: 8sig
#SL1: 2
#SL2: 2
#SL4: 0
#VB3R2: 49

#phyper((number of significant genes in this cross that overlap with other crosses),
#(number of genes significant in other crosses), 
#(number of genes not significant in other crosses), 
#(number of genes significant in this cross that don't overlap with any other cross)) 

#AR2 (HS1) 3 overlap, 66-(11-3), 152-58, 11-3
1-(phyper(3, 58, 94, 8)) #not sig 0.361

#AR4 (HS3) 2 overlap, 66-(4-2), 152-64, 4-2
1-(phyper(2, 64, 88, 2)) #sig at 0

#AR5 (HS4) 5 overlap, 66-(8-5), 152-63, 8-5
1-(phyper(5, 63, 89, 3)) #sig at 0

#SL1 (MS1) 1 overlap, 66-(2-1), 152-65, 2-1
1-(phyper(1, 65, 87, 1)) #sig at 0

#VB3R2 (LS3) 7 overlap, 66-(49-7), 152-24, 49-7
1-(phyper(7, 24, 128, 42)) #not sig at 0.325
