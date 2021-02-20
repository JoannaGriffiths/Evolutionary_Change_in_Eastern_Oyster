# https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
#https://www.rdocumentation.org/packages/qqman/versions/0.1.2/topics/manhattan

library("qqman")

setwd("~/Desktop/Oyster_exome")
setwd("~/OneDrive/Documents/exome_desktop/fet_sliding_window")
setwd("~/exome_desktop/")

Fst_means_AR2 = read.table("AR2-SE_sig_edit3.fet")
colnames(Fst_means_AR2) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe

Fst_means_AR3 = read.table("AR3-SE_sig_edit3.fet")
colnames(Fst_means_AR3) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe

Fst_means_AR4 = read.table("AR4-SE_sig_edit3.fet")
colnames(Fst_means_AR4) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe

Fst_means_AR5 = read.table("AR5-SE_sig_edit3.fet")
colnames(Fst_means_AR5) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe

Fst_means_VB2 = read.table("VB2-SE_sig_edit3.fet")
colnames(Fst_means_VB2) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe

Fst_means_VB3 = read.table("VB3-SE_sig_edit3.fet")
colnames(Fst_means_VB3) <- cbind("CHR","spp","gene","SNP","IDK1","IDK2","IDK3","P") #adding column names to dataframe


quartz()
manhattan(Fst_means_AR2, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)
manhattan(Fst_means_AR3, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)
manhattan(Fst_means_AR4, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)
manhattan(Fst_means_AR5, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)
manhattan(Fst_means_VB2, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)
manhattan(Fst_means_VB3, chr="CHR", bp="SNP", snp="IDK3", p="P", suggestiveline = 3.48, logp = F)

head(Fst_means$CHR)


setwd("~/Desktop/Oyster_exome/fet_sliding_window")

Fet_means_AR2 = read.table("VB3-SE_w20000_cov10-200_Manhattan.fet", header = T)

quartz()
manhattan(Fet_means_AR2, chr="CHR", bp="BP_middle", snp="SNP", p="fet", suggestiveline = F, genomewideline = F, logp = F, ylim =c(0,300))

Fet_means_AR2 = read.table("AR_VB_AVE-SE_w20000_cov10-200_Manhattan.fet", header = T)
quartz()
manhattan(Fet_means_AR2, chr="CHR", bp="BP_middle", snp="SNP", p="AR_ave", suggestiveline = F, genomewideline = F, logp = F, ylim =c(0,300))


setwd("~/Desktop/Oyster_exome")

Fet_means_AR2 = read.table("AR2-SE_final2.fet", header = T)

quartz()
manhattan(Fet_means_AR2, chr="CHR", bp="BP", snp="gene", p="fet", suggestiveline = F, genomewideline = F, logp = F, ylim =c(0,10))


setwd("~/Desktop/Oyster_exome/fet_sliding_window")
Results_VB3R2 = read.table("VB3R2-newpvalue_Manhattan", header = T)
windows()
manhattan(Results, chr="chr", bp="BP_position", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim=c(0,10), highlight=sig_list)

###start here
setwd("~/exome_desktop/fet_sliding_window")
Results_VB3R2 = read.table("VB3R2-newpvalue_Manhattan", header = T) #LA3
Results_VB3R1 = read.table("VB3R1-newpvalue_Manhattan", header = T) #LA2
Results_VB2 = read.table("VB2-newpvalue_Manhattan", header = T) #LA1
Results_AR2 = read.table("AR2-newpvalue_Manhattan", header = T) #TX1
Results_AR3 = read.table("AR3-newpvalue_Manhattan", header = T) #TX2
Results_AR4 = read.table("AR4-newpvalue_Manhattan", header = T) #TX3
Results_AR5 = read.table("AR5-newpvalue_Manhattan", header = T) #TX4
Results_SL1 = read.table("SL1-newpvalue_Manhattan", header = T) #LA4
Results_SL2 = read.table("SL2-newpvalue_Manhattan", header = T) #LA5
Results_SL3 = read.table("SL3-newpvalue_Manhattan", header = T) #LA6
Results_SL3NS = read.table("SL3-NS-newpvalue_Manhattan", header = T) #LA6 NS
#Results_SL4 = read.table("SL4-newpvalue_Manhattan", header = T)
#Results_SL4-NS = read.table("SL4-NS-newpvalue_Manhattan", header = T)

#sig_p = read.table("VB3R2-sigSNPs_forManhattan")
#getting allele freq
al_freq = read.delim("../SL1-SE_exact_rc_pop9_edit6", header=F)
al_freq<- al_freq[,c(2,19)]
colnames(al_freq) <- c('ID', 'freq')
Results$ID <- paste(Results$Gene_Pos_start, Results$Gene_Pos_end, Results$SNP, sep="-")
Results_alfreq <- merge(Results, al_freq, by='ID')
Results_alfreq$pos_freq <- (Results_alfreq$freq * -1) +1
#subset results to plot only sig allele freq changes
sig_Results_alfreq <- subset(Results_alfreq, padj < 0.05) #20

#AR2
AR2_sig <- subset(Results_AR2, padj < 0.05)
AR2_sig2 <- subset(AR2_sig, Gene=="22430406"| Gene=="22436694" | Gene=="22436855"| Gene=="22445237"| Gene=="22452549"| Gene=="22458334"| Gene=="22460042"| Gene=="22465298"| Gene=="22478344"| Gene=="22482159"| Gene=="22483747")
sig_list_AR2 <- AR2_sig2[,11]

#AR4
AR4_sig <- subset(Results_AR4, padj < 0.05)
AR4_sig2 <- subset(AR4_sig, Gene=="22466570"| Gene=="22469890" | Gene=="22477495"| Gene=="22482159")
sig_list_AR4 <- AR4_sig2[,11]

#AR5
AR5_sig <- subset(Results_AR5, padj < 0.05)
AR5_sig2 <- subset(AR5_sig, Gene=="22430406"| Gene=="22431296" | Gene=="22436490"| Gene=="22436855"| Gene=="22442587"| Gene=="22452045"| Gene=="22464577"| Gene=="22487619")
sig_list_AR5 <- AR5_sig2[,11]

#SL1
SL1_sig <- subset(Results_SL1, padj < 0.05)
SL1_sig2 <- subset(SL1_sig, Gene=="22447067"| Gene=="22480648")
sig_list_SL1 <- SL1_sig2[,11]

#SL2
SL2_sig <- subset(Results_SL2, padj < 0.05)
SL2_sig2 <- subset(SL2_sig, Gene=="22440183"| Gene=="22449816")
sig_list_SL2 <- SL2_sig2[,11]

#VB3_R2
VB3_R2_sig <- subset(Results_VB3R2, padj < 0.05)
VB3_R2_sig2 <- subset(VB3_R2_sig, gene=="22430406"| gene=="22430888"| gene=="22432300"| gene=="22434541"| gene=="22434694"| gene=="22434894"| gene=="22435092"| gene=="22436563"| gene=="22436855"| gene=="22439132"| gene=="22441293"| gene=="22442587"| gene=="22443236"| gene=="22446194"| gene=="22447067"| gene=="22448581"| gene=="22449029"| gene=="22450779"| gene=="22451612"| gene=="22451937"| gene=="22452045"| gene=="22452319"| gene=="22452520"| gene=="22455398"| gene=="22455778"| gene=="22455865"| gene=="22456272"| gene=="22456387"| gene=="22457334"| gene=="22457562"| gene=="22460685"| gene=="22462887"| gene=="22463410"| gene=="22464577"| gene=="22464888"| gene=="22464985"| gene=="22467005"| gene=="22468800"| gene=="22469890"| gene=="22470155"| gene=="22470793"| gene=="22473041"| gene=="22473715"| gene=="22480141"| gene=="22483345"| gene=="22484795"| gene=="22485000"| gene=="22486585"| gene=="22489168")
sig_list_VB3R2 <- VB3_R2_sig2[,12]

#plot all figures together, own axes
windows() #I think karyoploteR needs to be turned off for this to work? being masked by it
par(mfrow=c(3,2)) # 6 figures arranged in 3 rows and 2 columns
manhattan(Results_AR2, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_AR2)
manhattan(Results_AR4, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_AR4)
manhattan(Results_AR5, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_AR5)
manhattan(Results_SL1, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_SL1)
manhattan(Results_SL2, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_SL2)
manhattan(Results_VB3R2, chr="chr", bp="BP_position", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list_VB3R2)


#plot all figures with same x axis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
library("karyoploteR") #this didn't work when I wasn't connected to internet, really frustrating...

#help making GRanges object: https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/makeGRangesFromDataFrame.html
#https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html
AR2_GRange <- makeGRangesFromDataFrame(Results_AR2, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR2_GRange$pvalue <- Results_AR2$padj
AR2_GRange$Gene <- as.numeric(Results_AR2$Gene)
AR2_GRange

AR4_GRange <- makeGRangesFromDataFrame(Results_AR4, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR4_GRange$pvalue <- Results_AR4$padj
AR4_GRange$Gene <- as.numeric(Results_AR4$Gene)
AR4_GRange

AR5_GRange <- makeGRangesFromDataFrame(Results_AR5, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR5_GRange$pvalue <- Results_AR5$padj
AR5_GRange$Gene <- as.numeric(Results_AR5$Gene)
AR5_GRange

SL1_GRange <- makeGRangesFromDataFrame(Results_SL1, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL1_GRange$pvalue <- Results_SL1$padj
SL1_GRange$Gene <- as.numeric(Results_SL1$Gene)
SL1_GRange

SL2_GRange <- makeGRangesFromDataFrame(Results_SL2, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL2_GRange$pvalue <- Results_SL2$padj
SL2_GRange$Gene <- as.numeric(Results_SL2$Gene)
SL2_GRange

VB3R2_GRange <- makeGRangesFromDataFrame(Results_VB3R2, seqnames.field="chr", start.field="BP_position", end.field="BP_position")
VB3R2_GRange$pvalue <- Results_VB3R2$padj
VB3R2_GRange$Gene <- as.numeric(Results_VB3R2$gene)
VB3R2_GRange

SL3_GRange <- makeGRangesFromDataFrame(Results_SL3, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL3_GRange$pvalue <- Results_SL3$padj
SL3_GRange$Gene <- as.numeric(Results_SL3$Gene)

SL3NS_GRange <- makeGRangesFromDataFrame(Results_SL3NS, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
SL3NS_GRange$pvalue <- Results_SL3NS$padj
SL3NS_GRange$Gene <- as.numeric(Results_SL3NS$Gene)

AR3_GRange <- makeGRangesFromDataFrame(Results_AR3, seqnames.field="CHR", start.field="SNP_Pos", end.field="SNP_Pos")
AR3_GRange$pvalue <- Results_AR3$padj
AR3_GRange$Gene <- as.numeric(Results_AR3$Gene)
AR3_GRange

VB3R1_GRange <- makeGRangesFromDataFrame(Results_VB3R1, seqnames.field="CHR", start.field="SNP_pos", end.field="SNP_pos")
VB3R1_GRange$pvalue <- Results_VB3R1$padj
VB3R1_GRange$Gene <- as.numeric(Results_VB3R1$Gene)
VB3R1_GRange

VB2_GRange <- makeGRangesFromDataFrame(Results_VB2, seqnames.field="CHR", start.field="SNP_pos", end.field="SNP_pos")
VB2_GRange$pvalue <- Results_VB2$padj
VB2_GRange$Gene <- as.numeric(Results_VB2$Gene)
VB2_GRange

#create custom genome: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/CustomGenomes/CustomGenomes.html
c.virginica.genome <- toGRanges("cvirginica_exome_genome.txt")
#plotting help: https://rdrr.io/bioc/karyoploteR/man/kpPlotManhattan.html

AR2_sig <- subset(AR2_GRange, pvalue < 0.05)
AR2_sig_snps <- subset(AR2_sig, Gene=="22430406"| Gene=="22436694" | Gene=="22436855"| Gene=="22445237"| Gene=="22452549"| Gene=="22458334"| Gene=="22460042"| Gene=="22465298"| Gene=="22478344"| Gene=="22482159"| Gene=="22483747")
AR2_sig_snps

AR4_sig <- subset(AR4_GRange, pvalue < 0.05)
AR4_sig_snps <- subset(AR4_sig, Gene=="22466570"| Gene=="22469890" | Gene=="22477495"| Gene=="22482159")
AR4_sig_snps

AR5_sig <- subset(AR5_GRange, pvalue < 0.05)
AR5_sig_snps <- subset(AR5_sig, Gene=="22430406"| Gene=="22431296" | Gene=="22436490"| Gene=="22436855"| Gene=="22442587"| Gene=="22452045"| Gene=="22464577"| Gene=="22487619")
AR5_sig_snps

SL1_sig <- subset(SL1_GRange, pvalue < 0.05)
SL1_sig_snps <- subset(SL1_sig, Gene=="22447067"| Gene=="22480648")
SL1_sig_snps

SL2_sig <- subset(SL2_GRange, pvalue < 0.05)
SL2_sig_snps <- subset(SL2_sig, Gene=="22440183"| Gene=="22449816")
SL2_sig_snps

VB3R2_sig <- subset(VB3R2_GRange, pvalue < 0.05)
VB3R2_sig_snps <- subset(VB3R2_sig, Gene=="22430406"| Gene=="22430888"| Gene=="22432300"| Gene=="22434541"| Gene=="22434694"| Gene=="22434894"| Gene=="22435092"| Gene=="22436563"| Gene=="22436855"| Gene=="22439132"| Gene=="22441293"| Gene=="22442587"| Gene=="22443236"| Gene=="22446194"| Gene=="22447067"| Gene=="22448581"| Gene=="22449029"| Gene=="22450779"| Gene=="22451612"| Gene=="22451937"| Gene=="22452045"| Gene=="22452319"| Gene=="22452520"| Gene=="22455398"| Gene=="22455778"| Gene=="22455865"| Gene=="22456272"| Gene=="22456387"| Gene=="22457334"| Gene=="22457562"| Gene=="22460685"| Gene=="22462887"| Gene=="22463410"| Gene=="22464577"| Gene=="22464888"| Gene=="22464985"| Gene=="22467005"| Gene=="22468800"| Gene=="22469890"| Gene=="22470155"| Gene=="22470793"| Gene=="22473041"| Gene=="22473715"| Gene=="22480141"| Gene=="22483345"| Gene=="22484795"| Gene=="22485000"| Gene=="22486585"| Gene=="22489168")
VB3R2_sig_snps

overlap1 <- subset(AR5_sig_snps, Gene=="22464577"| Gene=="22430406"|Gene=="22436855" |Gene=="22442587"|Gene=="22452045")
overlap2 <- subset(AR4_sig_snps, Gene=="22482159"|Gene=="22469890")
overlap3 <- subset(VB3R2_sig_snps, Gene=="22447067")
reg1 <- extendRegions(overlap1, 30e5, 30e5)
reg2 <- extendRegions(overlap2, 30e5, 30e5)
reg3 <- extendRegions(overlap3, 30e5, 30e5)

##Manhattan plots for non-significant crosses
grDevices::windows()
kp <- plotKaryotype(genome = c.virginica.genome, plot.type = 4, cex=0.8)
kp <- kpPlotManhattan(kp, data=VB2_GRange, pval = VB2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(5,5),
                      genomewideline =0,
                      ymin=0, ymax=8
)

kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(5,5), cex=0.5)
kpAddLabels(kp, labels = "LA1", side="right", r0=autotrack(5,5), cex=0.8)
kp <- kpPlotManhattan(kp, data=VB3R1_GRange, pval = VB3R1_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(4,5),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(4,5), cex=0.5)
kpAddLabels(kp, labels = "LA2", side="right", r0=autotrack(4,5), cex=0.8)
kp <- kpPlotManhattan(kp, data=SL3_GRange, pval = SL3_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(3,5),
                      genomewideline =0,
                      ymin=0, ymax=8
                      )
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(3,5), cex=0.5)
kpAddLabels(kp, labels = "LA6", side="right", r0=autotrack(3,5), cex=0.8)
kp <- kpPlotManhattan(kp, data=SL3NS_GRange, pval = SL3NS_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(2,5),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(2,5), cex=0.5)
kpAddLabels(kp, labels = "LA6 NS",side="right", r0=autotrack(2,5), cex=0.8)
kp <- kpPlotManhattan(kp, data=AR3_GRange, pval = AR3_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      r0=autotrack(1,5),
                      genomewideline =0,
                      ymin=0, ymax=8
)
kpAxis(kp, ymin=0, ymax=8, tick.pos = c(2,4,6), r0=autotrack(1,5), cex=0.5)
kpAddLabels(kp, labels = "TX2", side="right", r0=autotrack(1,5), cex=0.8)
kpAddLabels(kp, labels = "-log(pvalue)", srt=90, label.margin = 0.04)




##Manhattan plots for significant crosses
grDevices::windows()
kp <- plotKaryotype(genome = c.virginica.genome, plot.type = 4, cex=0.8)
kp <- kpPlotManhattan(kp, data=AR2_GRange, pval = AR2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR2_sig_snps,
                      r0=autotrack(3,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(3,6), cex=0.5)
kpAddLabels(kp, labels = "TX1", side="right", r0=autotrack(3,6), cex=0.8)
kp <- kpPlotManhattan(kp, data=AR4_GRange, pval = AR4_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR4_sig_snps,
                      r0=autotrack(2,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(2,6), cex=0.5)
kpAddLabels(kp, labels = "TX3",side="right", r0=autotrack(2,6), cex=0.8)
kp <- kpPlotManhattan(kp, data=AR5_GRange, pval = AR5_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = AR5_sig_snps,
                      r0=autotrack(1,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(1,6), cex=0.5)
kpAddLabels(kp, labels = "TX4", side="right", r0=autotrack(1,6), cex=0.8)
kp <- kpPlotManhattan(kp, data=SL1_GRange, pval = SL1_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = SL1_sig_snps,
                      r0=autotrack(5,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(5,6), cex=0.5)
kpAddLabels(kp, labels = "LA4", side="right", r0=autotrack(5,6), cex=0.8)
kp <- kpPlotManhattan(kp, data=SL2_GRange, pval = SL2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = SL2_sig_snps,
                      r0=autotrack(4,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(4,6), cex=0.5)
kpAddLabels(kp, labels = "LA5", side="right", r0=autotrack(4,6), cex=0.8)
kp <- kpPlotManhattan(kp, data=VB3R2_GRange, pval = VB3R2_GRange$pvalue, 
                      suggestive.col="blue", suggestive.lwd = 1.3, suggestiveline = 1.3, 
                      points.cex = 0.5, 
                      highlight = VB3R2_sig_snps,
                      r0=autotrack(6,6),
                      genomewideline =0,
                      ymin=0, ymax=10
)
kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(6,6), cex=0.5)
kpAddLabels(kp, labels = "LA3", side="right", r0=autotrack(6,6), cex=0.8)
kpAddLabels(kp, labels = "-log(pvalue)", srt=90, label.margin = 0.04)
kpRect(kp, data=reg1, y0=0, y1=1, col=NA, border="red", lwd=1)
kpRect(kp, data=reg2, y0=0, y1=1, col=NA, border="red", lwd=1)
kpRect(kp, data=reg3, y0=0, y1=1, col=NA, border="red", lwd=1)

#kpAxis(kp, ymin=0, ymax=10, tick.pos = c(2,4,6,8), r0=autotrack(3,6), cex=0.5)
#kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax, tick.pos = c(0, seq_len(floor(kp$latest.plot$computed.values$ymax))), r0=autotrack(2,6), cex=0.5)


#plot manhattan and allele freq on same graph
windows()
par(mfrow=c(3,2))
manhattan(Results_AR2, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list)
par(new=T)
manhattan(sig_Results_alfreq, chr="CHR", bp="SNP_Pos", snp="SNP_Pos", p="pos_freq", suggestiveline = F, genomewideline = F, logp = F, type = "p", axes = FALSE, bty = "n", xlab = "", ylab = "", col = 'red')
axis(side=4, at = pretty(range(1,1.3)))
mtext("z", side=4, line=3)

####test to figure out locations of significant genes overlapping between crosses
sig <- subset(Results, padj < 0.05)
sig2 <- subset(sig, Gene=="22430406")
sig_list <- sig2[,11]

manhattan(Results, chr="CHR", bp="SNP_Pos", snp="Man_SNP", p="padj", suggestiveline = 1.3, genomewideline = F, logp = T, ylim =c(0,10), highlight = sig_list)

