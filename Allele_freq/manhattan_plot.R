

library("dplyr")
setwd("~/LSU/Research/Oyster Exome Capture Experiment")

### read in annotation file for each gene
annot_final2 <- read.delim("Supplemental_Table1", header = T)

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



###################
#plot all figures with same x axis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
library("karyoploteR") #this didn't work when I wasn't connected to internet, really frustrating...

#help making GRanges object: https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/makeGRangesFromDataFrame.html
#https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html
AR2_GRange <- makeGRangesFromDataFrame(AR2, seqnames.field="CHR", start.field="SNP_pos", end.field="SNP_pos")
AR2_GRange$pvalue <- AR2$padj
AR2_GRange$Gene <- as.numeric(AR2$Gene)
AR2_GRange

AR4_GRange <- makeGRangesFromDataFrame(AR4, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR4_GRange$pvalue <- AR4$padj
AR4_GRange$Gene <- as.numeric(AR4$Gene)
AR4_GRange

AR5_GRange <- makeGRangesFromDataFrame(AR5, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR5_GRange$pvalue <- AR5$padj
AR5_GRange$Gene <- as.numeric(AR5$Gene)
AR5_GRange

SL1_GRange <- makeGRangesFromDataFrame(SL1, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL1_GRange$pvalue <- SL1$padj
SL1_GRange$Gene <- as.numeric(SL1$Gene)
SL1_GRange

VB2_GRange <- makeGRangesFromDataFrame(VB2, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
VB2_GRange$pvalue <- VB2$padj
VB2_GRange$Gene <- as.numeric(VB2$Gene)
VB2_GRange

VB3R2_GRange <- makeGRangesFromDataFrame(VB3_R2, seqnames.field="chr", start.field="SNP_pos", end.field="SNP_pos")
VB3R2_GRange$pvalue <- VB3_R2$padj
VB3R2_GRange$Gene <- as.numeric(VB3_R2$Gene)
VB3R2_GRange

SL3_GRange <- makeGRangesFromDataFrame(SL3, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL3_GRange$pvalue <- SL3$padj
SL3_GRange$Gene <- as.numeric(SL3$Gene)

SL3NS_GRange <- makeGRangesFromDataFrame(SLNS3, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL3NS_GRange$pvalue <- SLNS3$padj
SL3NS_GRange$Gene <- as.numeric(SLNS3$Gene)

AR3_GRange <- makeGRangesFromDataFrame(AR3, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR3_GRange$pvalue <- AR3$padj
AR3_GRange$Gene <- as.numeric(AR3$Gene)
AR3_GRange

VB3R1_GRange <- makeGRangesFromDataFrame(VB3, seqnames.field="CHR", start.field="SNP_pos", end.field="SNP_pos")
VB3R1_GRange$pvalue <- VB3$padj
VB3R1_GRange$Gene <- as.numeric(VB3$Gene)
VB3R1_GRange

SL2_GRange <- makeGRangesFromDataFrame(SL2, seqnames.field="CHR", start.field="SNP_pos", end.field="SNP_pos")
SL2_GRange$pvalue <- SL2$padj
SL2_GRange$Gene <- as.numeric(SL2$Gene)
SL2_GRange

#create custom genome: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/CustomGenomes/CustomGenomes.html
c.virginica.genome <- toGRanges("~/exome_desktop/fet_sliding_window/cvirginica_exome_genome.txt")
#plotting help: https://rdrr.io/bioc/karyoploteR/man/kpPlotManhattan.html

AR2_sig <- subset(AR2_GRange, pvalue < 0.05)
AR2_sig_snps <- subset(AR2_sig, Gene=="22436855"| Gene=="22452549" | Gene=="22458334"| Gene=="22476633")
AR2_sig_snps

AR4_sig <- subset(AR4_GRange, pvalue < 0.05)
AR4_sig_snps <- subset(AR4_sig, Gene=="22446194"| Gene=="22460685" | Gene=="22464985"| Gene=="22473715")
AR4_sig_snps

AR5_sig <- subset(AR5_GRange, pvalue < 0.05)
AR5_sig_snps <- subset(AR5_sig, Gene=="22430406"| Gene=="22430468" | Gene=="22430888"| Gene=="22436272"| Gene=="22439003"| Gene=="22447067"| Gene=="22450779"| Gene=="22455865"| Gene=="22456387"| Gene=="22457562"| Gene=="22458334"| Gene=="22464985"| Gene=="22477495")
AR5_sig_snps

SL1_sig <- subset(SL1_GRange, pvalue < 0.05)
SL1_sig_snps <- subset(SL1_sig, Gene=="22447067"| Gene=="22477509")
SL1_sig_snps

VB2_sig <- subset(VB2_GRange, pvalue < 0.05)
VB2_sig_snps <- subset(VB2_sig, Gene=="22441408"| Gene=="22460202" | Gene=="22467005")
VB2_sig_snps

VB3R1_sig <- subset(VB3R1_GRange, pvalue < 0.05)
VB3R1_sig_snps <- subset(VB3R1_sig, Gene=="22460042")
VB3R1_sig_snps

VB3R2_sig <- subset(VB3R2_GRange, pvalue < 0.05)
VB3R2_sig_snps <- subset(VB3R2_sig, Gene=="22430406"| Gene=="22430888"| Gene=="22434940"| Gene=="22436694"| Gene=="22436855"| Gene=="22437155"| Gene=="22439132"| Gene=="22443236"| Gene=="22446194"| Gene=="22447067"| Gene=="22449029"| Gene=="22452045"| Gene=="22452419"| Gene=="22452520"| Gene=="22454425"| Gene=="22455778"| Gene=="22455865"| Gene=="22456272"| Gene=="22458334"| Gene=="22464888"| Gene=="22464985"| Gene=="22465298"| Gene=="22467005"| Gene=="22473715"| Gene=="22475843"| Gene=="22477495"| Gene=="22479437"| Gene=="22484795"| Gene=="22486585"| Gene=="22486637")
VB3R2_sig_snps

overlap1 <- subset(AR5_sig_snps, Gene=="22430406"| Gene=="22430888"|Gene=="22447067" |Gene=="22455865"|Gene=="22458334")
overlap2 <- subset(VB3R2_sig_snps, Gene=="22436855" |Gene=="22446194"|Gene=="22467005"|Gene=="22477495")
overlap3 <- subset(AR5_sig_snps, Gene=="22464985")
overlap4 <- subset(AR4_sig_snps, Gene=="22473715")
reg1 <- extendRegions(overlap1, 30e5, 30e5)
reg2 <- extendRegions(overlap2, 30e5, 30e5)
reg3 <- extendRegions(overlap3, 50e5, 50e5)
#reg4 <- extendRegions(overlap4, 50e5, 50e5)


##Manhattan plots for non-significant crosses
grDevices::windows()
kp <- plotKaryotype(genome = c.virginica.genome, plot.type = 4, cex=0.8)

kp <- kpPlotManhattan(kp, data=SL2_GRange, pval = SL2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(4,4),
                      genomewideline =0,
                      ymin=0, ymax=8
)

kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(4,4), cex=0.5)
kpAddLabels(kp, labels = "LA5", side="right", r0=autotrack(4,4), cex=0.8)

kp <- kpPlotManhattan(kp, data=SL3_GRange, pval = SL3_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(3,4),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(3,4), cex=0.5)
kpAddLabels(kp, labels = "LA6", side="right", r0=autotrack(3,4), cex=0.8)

kp <- kpPlotManhattan(kp, data=SL3NS_GRange, pval = SL3NS_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(2,4),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(2,4), cex=0.5)
kpAddLabels(kp, labels = "LA6 NS",side="right", r0=autotrack(2,4), cex=0.8)

kp <- kpPlotManhattan(kp, data=AR3_GRange, pval = AR3_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(1,4),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(1,4), cex=0.5)
kpAddLabels(kp, labels = "TX2", side="right", r0=autotrack(1,4), cex=0.8)
kpAddLabels(kp, labels = "-log(pvalue)", srt=90, label.margin = 0.04)




##Manhattan plots for significant crosses
grDevices::windows()
kp <- plotKaryotype(genome = c.virginica.genome, plot.type = 4, cex=0.8)

kp <- kpPlotManhattan(kp, data=VB3R1_GRange, pval = VB3R1_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = VB3R1_sig_snps,
                      r0=autotrack(6,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(6,7), cex=0.5)
kpAddLabels(kp, labels = "LA2", side="right", r0=autotrack(6,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=AR2_GRange, pval = AR2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR2_sig_snps,
                      r0=autotrack(3,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(3,7), cex=0.5)
kpAddLabels(kp, labels = "TX1", side="right", r0=autotrack(3,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=AR4_GRange, pval = AR4_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR4_sig_snps,
                      r0=autotrack(2,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(2,7), cex=0.5)
kpAddLabels(kp, labels = "TX3",side="right", r0=autotrack(2,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=AR5_GRange, pval = AR5_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR5_sig_snps,
                      r0=autotrack(1,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(1,7), cex=0.5)
kpAddLabels(kp, labels = "TX4", side="right", r0=autotrack(1,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=SL1_GRange, pval = SL1_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = SL1_sig_snps,
                      r0=autotrack(4,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(4,7), cex=0.5)
kpAddLabels(kp, labels = "LA4", side="right", r0=autotrack(4,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=VB2_GRange, pval = VB2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = VB2_sig_snps,
                      r0=autotrack(7,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(7,7), cex=0.5)
kpAddLabels(kp, labels = "LA1", side="right", r0=autotrack(7,7), cex=0.8)

kp <- kpPlotManhattan(kp, data=VB3R2_GRange, pval = VB3R2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = VB3R2_sig_snps,
                      r0=autotrack(5,7),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(5,7), cex=0.5)
kpAddLabels(kp, labels = "LA3", side="right", r0=autotrack(5,7), cex=0.8)
kpAddLabels(kp, labels = "-log(pvalue)", srt=90, label.margin = 0.04)
kpRect(kp, data=reg1, y0=0, y1=1, col=NA, border="red", lwd=1)
kpRect(kp, data=reg2, y0=0, y1=1, col=NA, border="red", lwd=1)
kpRect(kp, data=reg3, y0=0, y1=1, col=NA, border="red", lwd=1)


