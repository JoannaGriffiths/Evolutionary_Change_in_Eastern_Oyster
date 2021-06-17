


setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("AR2-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 15)]

newpvalue = read.delim("AR2-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #92
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #74
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #9
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #9

sig_gene_overlap <- subset(sig_newpvalue_overlap, Orientation == "overlapping")
write.table(sig_gene_overlap, file = "AR2-sig_gene_overlap.txt", sep = "\t", quote=F)

setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("AR3-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 16)]

newpvalue = read.delim("AR3-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

newpvalue_overlap$padj <- as.numeric(as.character(newpvalue_overlap$padj))
sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #10
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #6
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #3
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #1


setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("AR4-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 17)]

newpvalue = read.delim("AR4-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

newpvalue_overlap$padj <- as.numeric(as.character(newpvalue_overlap$padj))
sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #60
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #41
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #6
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #13



setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("AR5-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 18)]

newpvalue = read.delim("AR5-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

newpvalue_overlap$padj <- as.numeric(as.character(newpvalue_overlap$padj))
sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #103
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #67
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #20
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #16



setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("SL1-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 19)]

newpvalue = read.delim("SL1-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

newpvalue_overlap$padj <- as.numeric(as.character(newpvalue_overlap$padj))
sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #24
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #13
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #6
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #5


setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("SL2-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 19)]

newpvalue = read.delim("SL2-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #1
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #0
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #1
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #0


setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("VB3R1-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 16)]

newpvalue = read.delim("VB3R1-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$CHR, newpvalue$SNP_Pos, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #20
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #12
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #5
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #3




setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Overlap = read.delim("VB3R2-geneoverlap.txt", header = T)
Overlap$Merge_ID <- paste(Overlap$Chromosome, Overlap$Start, sep="-")
Overlap2 <- Overlap[,c(12, 13, 17)]

newpvalue = read.delim("VB3R2-newpvalue_SeqMonk", header = T)
newpvalue$Merge_ID <- paste(newpvalue$chr, newpvalue$BP_position, sep="-")

newpvalue_overlap <- merge(newpvalue, Overlap2, by='Merge_ID')

newpvalue_overlap$padj <- as.numeric(as.character(newpvalue_overlap$padj))
sig_newpvalue_overlap <- subset(newpvalue_overlap, padj < 0.05) #339
length(which(sig_newpvalue_overlap$Orientation == "overlapping")) #225
length(which(sig_newpvalue_overlap$Orientation == "upstream")) #70
length(which(sig_newpvalue_overlap$Orientation == "downstream")) #44


setwd("C:/Users/joann/Desktop/oyster_exome")
AR2_SNPs <- read.delim("AR2-start_SNPs.txt", header = T)
length(which(AR2_SNPs$Orientation == "overlapping")) #4219
length(which(AR2_SNPs$Orientation == "upstream")) #692
length(which(AR2_SNPs$Orientation == "downstream")) #816

AR3_SNPs <- read.delim("AR3-start_SNPs.txt", header = T)
length(which(AR3_SNPs$Orientation == "overlapping")) #4281
length(which(AR3_SNPs$Orientation == "upstream")) #1345
length(which(AR3_SNPs$Orientation == "downstream")) #1234

AR4_SNPs <- read.delim("AR4-start_SNPs.txt", header = T)
length(which(AR4_SNPs$Orientation == "overlapping")) #5704
length(which(AR4_SNPs$Orientation == "upstream")) #1196
length(which(AR4_SNPs$Orientation == "downstream")) #1275

AR5_SNPs <- read.delim("AR5-start_SNPs.txt", header = T)
length(which(AR5_SNPs$Orientation == "overlapping")) #7877
length(which(AR5_SNPs$Orientation == "upstream")) #1794
length(which(AR5_SNPs$Orientation == "downstream")) #1724

SL1_SNPs <- read.delim("SL1-start_SNPs.txt", header = T)
length(which(SL1_SNPs$Orientation == "overlapping")) #3474
length(which(SL1_SNPs$Orientation == "upstream")) #843
length(which(SL1_SNPs$Orientation == "downstream")) #943

SL2_SNPs <- read.delim("SL2-start_SNPs.txt", header = T)
length(which(SL2_SNPs$Orientation == "overlapping")) #3206
length(which(SL2_SNPs$Orientation == "upstream")) #820
length(which(SL2_SNPs$Orientation == "downstream")) #823

VB2_SNPs <- read.delim("VB2-start_SNPs.txt", header = T)
length(which(VB2_SNPs$Orientation == "overlapping")) #5749
length(which(VB2_SNPs$Orientation == "upstream")) #1475
length(which(VB2_SNPs$Orientation == "downstream")) #1484

VB3R1_SNPs <- read.delim("VB3-R1-start_SNPs.txt", header = T)
length(which(VB3R1_SNPs$Orientation == "overlapping")) #5722
length(which(VB3R1_SNPs$Orientation == "upstream")) #1193
length(which(VB3R1_SNPs$Orientation == "downstream")) #1149

VB3R2_SNPs <- read.delim("VB3-R2-start_SNPs.txt", header = T)
length(which(VB3R2_SNPs$Orientation == "overlapping")) #4219
length(which(VB3R2_SNPs$Orientation == "upstream")) #692
length(which(VB3R2_SNPs$Orientation == "downstream")) #816

AR2_df <- data.frame(
  start=c(74,9,9),
  end=c(4219,682,816))
chisq.test(AR2_df) #not sig 0.33

AR3_df <- data.frame(
  start=c(6,3,1),
  end=c(4281,1345,1234))
chisq.test(AR3_df) #not sig 0.63

AR4_df <- data.frame(
  start=c(41,6,13),
  end=c(5704,1196,1275))
chisq.test(AR4_df) #not sig 0.32

AR5_df <- data.frame(
  start=c(67,20,16),
  end=c(7877,1794,1724))
chisq.test(AR5_df) #not sig 0.56

SL1_df <- data.frame(
  start=c(13,6,5),
  end=c(3474,843,943))
chisq.test(SL1_df) #not sig 0.40

SL2_df <- data.frame(
  start=c(0,1,0),
  end=c(3206,820,823))
chisq.test(SL2_df) #not sig 0.086

VB2_df <- data.frame(
  start=c(0,1,0),
  end=c(5749,1475,1484))
chisq.test(VB2_df) #not sig 0.086

