


setwd("~/LSU/Research/Oyster Exome Capture Experiment/fet/June18_2021")
AR2_start_Overlap = read.delim("Annotated Probe Report for AR2_start.txt", header = T)
AR2_sig_Overlap = read.delim("Annotated Probe Report for AR2_sig_SNPs.txt", header = T)

length(which(AR2_start_Overlap$Feature.Orientation == "overlapping")) #2016
length(which(AR2_start_Overlap$Feature.Orientation == "upstream")) #299
length(which(AR2_start_Overlap$Feature.Orientation == "downstream")) #357

length(which(AR2_sig_Overlap$Feature.Orientation == "overlapping")) #30
length(which(AR2_sig_Overlap$Feature.Orientation == "upstream")) #5
length(which(AR2_sig_Overlap$Feature.Orientation == "downstream")) #10

AR2_df <- data.frame(
  start=c(2016,299,357),
  sig=c(30,5,10))
chisq.test(AR2_df) #not sig 0.22


AR4_start_Overlap = read.delim("Annotated Probe Report for AR4_start.txt", header = T)
AR4_sig_Overlap = read.delim("Annotated Probe Report for AR4_sig_SNPs.txt", header = T)

length(which(AR4_start_Overlap$Feature.Orientation == "overlapping")) #2714
length(which(AR4_start_Overlap$Feature.Orientation == "upstream")) #572
length(which(AR4_start_Overlap$Feature.Orientation == "downstream")) #527

length(which(AR4_sig_Overlap$Feature.Orientation == "overlapping")) #31
length(which(AR4_sig_Overlap$Feature.Orientation == "upstream")) #5
length(which(AR4_sig_Overlap$Feature.Orientation == "downstream")) #3

AR4_df <- data.frame(
  start=c(2714,572,527),
  sig=c(31,5,3))
chisq.test(AR4_df) #not sig 0.46


AR5_start_Overlap = read.delim("Annotated Probe Report for AR5_start.txt", header = T)
AR5_sig_Overlap = read.delim("Annotated Probe Report for AR5_sig_SNPs.txt", header = T)

length(which(AR5_start_Overlap$Feature.Orientation == "overlapping")) #3896
length(which(AR5_start_Overlap$Feature.Orientation == "upstream")) #925
length(which(AR5_start_Overlap$Feature.Orientation == "downstream")) #802

length(which(AR5_sig_Overlap$Feature.Orientation == "overlapping")) #79
length(which(AR5_sig_Overlap$Feature.Orientation == "upstream")) #17
length(which(AR5_sig_Overlap$Feature.Orientation == "downstream")) #19

AR5_df <- data.frame(
  start=c(3896,925,802),
  sig=c(79,17,19))
chisq.test(AR5_df) #not sig 0.74


SL1_start_Overlap = read.delim("Annotated Probe Report for SL1_start.txt", header = T)
SL1_sig_Overlap = read.delim("Annotated Probe Report for SL1_sig_SNPs.txt", header = T)

length(which(SL1_start_Overlap$Feature.Orientation == "overlapping")) #1456
length(which(SL1_start_Overlap$Feature.Orientation == "upstream")) #360
length(which(SL1_start_Overlap$Feature.Orientation == "downstream")) #354

length(which(SL1_sig_Overlap$Feature.Orientation == "overlapping")) #8
length(which(SL1_sig_Overlap$Feature.Orientation == "upstream")) #5
length(which(SL1_sig_Overlap$Feature.Orientation == "downstream")) #10

SL1_df <- data.frame(
  start=c(1456,360,354), #0.671, 0.1658, 0.1631
  sig=c(8,5,10)) #0.3478, 0.217, 0.4347
chisq.test(SL1_df) #sig 0.0008936


VB2_start_Overlap = read.delim("Annotated Probe Report for VB2_start.txt", header = T)
VB2_sig_Overlap = read.delim("Annotated Probe Report for VB2_sig_SNPs.txt", header = T)

length(which(VB2_start_Overlap$Feature.Orientation == "overlapping")) #2947
length(which(VB2_start_Overlap$Feature.Orientation == "upstream")) #701
length(which(VB2_start_Overlap$Feature.Orientation == "downstream")) #665

length(which(VB2_sig_Overlap$Feature.Orientation == "overlapping")) #32
length(which(VB2_sig_Overlap$Feature.Orientation == "upstream")) #3
length(which(VB2_sig_Overlap$Feature.Orientation == "downstream")) #1

VB2_df <- data.frame(
  start=c(2947,701,665), #0.6832, 0.1625, 0.1542
  sig=c(32,3,1)) #0.88, 0.0833, 0.0277 more likely to be found within gene body
chisq.test(VB2_df) #sig 0.0258


VB3_start_Overlap = read.delim("Annotated Probe Report for VB3_start.txt", header = T)
VB3_sig_Overlap = read.delim("Annotated Probe Report for VB3_sig_SNPs.txt", header = T)

length(which(VB3_start_Overlap$Feature.Orientation == "overlapping")) #2736
length(which(VB3_start_Overlap$Feature.Orientation == "upstream")) #540
length(which(VB3_start_Overlap$Feature.Orientation == "downstream")) #472

length(which(VB3_sig_Overlap$Feature.Orientation == "overlapping")) #22
length(which(VB3_sig_Overlap$Feature.Orientation == "upstream")) #5
length(which(VB3_sig_Overlap$Feature.Orientation == "downstream")) #4

VB3_df <- data.frame(
  start=c(2736,540,472),
  sig=c(22,5,4))
chisq.test(VB3_df) #not sig 0.95


VB3_R2_start_Overlap = read.delim("Annotated Probe Report for VB3_R2_start.txt", header = T)
VB3_R2_sig_Overlap = read.delim("Annotated Probe Report for VB3_R2_sig_SNPs.txt", header = T)

length(which(VB3_R2_start_Overlap$Feature.Orientation == "overlapping")) #1960
length(which(VB3_R2_start_Overlap$Feature.Orientation == "upstream")) #511
length(which(VB3_R2_start_Overlap$Feature.Orientation == "downstream")) #437

length(which(VB3_R2_sig_Overlap$Feature.Orientation == "overlapping")) #176
length(which(VB3_R2_sig_Overlap$Feature.Orientation == "upstream")) #31
length(which(VB3_R2_sig_Overlap$Feature.Orientation == "downstream")) #32

VB3_R2_df <- data.frame(
  start=c(1960,511,437),
  sig=c(176,31,32))
chisq.test(VB3_R2_df) #not sig 0.11
