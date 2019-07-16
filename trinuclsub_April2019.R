#NEDS

#DESEQ2-esque analysis
#each signature = a gene, each cluster = a condition

#1. Differential expression tables
#2. MA plot w/ highlights
#3. scatter for up and down regulated within one cluster + volcan for between two
#3. Boxplot + stats (anova etc?)
#4. Density plots (looks like porbablas)
#5. interaction plots (clusters with lines)

library(MutationalPatterns)
library(SomaticSignatures)
library(GenomicRanges)
library(ggplot2)
library(SomaticCancerAlterations)
library(NMF)
library(DESeq2)
library(edgeR)
library(pasilla)
library(dplyr)
library(apeglm)
#somaticsignatures decomposition methods = nmf (what we used before) or pca. we're aslso looking at trying a "custom" method - the nsNMF nonsmooth nonneg mat factorization


setwd("/Users/belfordak/Desktop/ML_mutation_viral_taxonomy/April_2019/")

#read in data
ML_table <- read.csv("./enrich_deplt_scores_VP1s.csv.csv", header = TRUE, row.names = 1, check.names = FALSE) #note, doesn't actually have clusters anymore
#merge same column names (clusts)
#ml_test <- t(rowsum(t(ML_table), group = colnames(ML_table), na.rm = T))
ml_test <- t(ML_table)
ml_test[is.na(ml_test)] <- 0
mlt <- ml_test * 1000
ml_test <- round(mlt)

#Need to:
#1. test of mutationalpatterns package
#2. test deseq *
#3. test edgeR

coldata = read.csv("./coldata_VP1s.csv", header = TRUE, row.names = 1, check.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = ml_test, colData = coldata , design = ~condition)
dds <- DESeq(dds)

res <- results(dds)
summary(res)
res

#actually producing all the data - pretty much have to do it like this
#clust 0 vs all
res_0_vs_1 <- results(dds, c('condition', "zero", 'one'))
res_0_vs_2 <- results(dds, c('condition', "zero", 'two'))
#res_0_vs_3 <- results(dds, c('condition', "zero", 'three'))
#res_0_vs_4 <- results(dds, c('condition', "zero", 'four'))
#res_0_vs_5 <- results(dds, c('condition', "zero", 'five'))
#clust 1 vs all
res_1_vs_2 <- results(dds, c('condition', "one", 'two'))
#res_1_vs_3 <- results(dds, c('condition', "one", 'three'))
#res_1_vs_4 <- results(dds, c('condition', "one", 'four'))
#res_1_vs_5 <- results(dds, c('condition', "one", 'five'))
#clust 2 vs all
#res_2_vs_3 <- results(dds, c('condition', "two", 'three'))
#res_2_vs_4 <- results(dds, c('condition', "two", 'four'))
#res_2_vs_5 <- results(dds, c('condition', "two", 'five'))
#clust 3 vs all
#res_3_vs_4 <- results(dds, c('condition', "three", 'four'))
#res_3_vs_5 <- results(dds, c('condition', "three", 'five'))
#clust 4 vs all
#res_4_vs_5 <- results(dds, c('condition', "four", 'five'))
#5 = nothing left, all comparisons done (just remember up/down regulartion inverse)

#out all to csv for excel manipulation
write.csv(res_0_vs_1, "./R_out/VP1s_res_0_vs_1.csv")
write.csv(res_0_vs_2, "./R_out/VP1s_res_0_vs_2.csv")
#write.csv(res_0_vs_3, "./R_out/FG_res_0_vs_3.csv")
#write.csv(res_0_vs_4, "./R_out/FG_res_0_vs_4.csv")
#write.csv(res_0_vs_5, "./R_out/res_0_vs_5.csv")

write.csv(res_1_vs_2, "./R_out/VP1s_res_1_vs_2.csv")
#write.csv(res_1_vs_3, "./R_out/FG_res_1_vs_3.csv")
#write.csv(res_1_vs_4, "./R_out/res_1_vs_4.csv")
#write.csv(res_1_vs_5, "./R_out/res_1_vs_5.csv")

write.csv(res_2_vs_3, "./R_out/FG_res_2_vs_3.csv")
#write.csv(res_2_vs_4, "./R_out/res_2_vs_4.csv")
#write.csv(res_2_vs_5, "./R_out/res_2_vs_5.csv")

#write.csv(res_3_vs_4, "./R_out/res_3_vs_4.csv")
#write.csv(res_3_vs_5, "./R_out/res_3_vs_5.csv")

#write.csv(res_4_vs_5, "./R_out/res_4_vs_5.csv")


#further DESEQ
plotMA(res_0_vs_5)


#MutationalPatterns
library(gridExtra)
#mlt <- t(ml_test)


#simple barplots of greatest logFC2 for each cluster comparison - I pulled this from the excel sheet
library(ggplot2)
#FG_0vs1
log_FC_2 = c(0.835098645,
             0.743653997,
             0.477778029,
             0.448687069,
             0.360648763,
             0.344085675,
             0.335484697,
             0.324464229,
             0.311930531,
             0.299650784)
sig_name = c("C[T>G]A",
             "C[C>G]A",
             "C[C>G]C",
             "C[T>G]C",
             "C[T>G]T",
             "G[T>A]C",
             "C[C>G]T",
             "T[C>G]G",
             "T[C>A]G",
             "G[C>T]G")

FG_0v1 <- data.frame(sig_name,log_FC_2)

p <- ggplot(data=FG_0v1, aes(x=sig_name, y=log_FC_2, fill = sig_name)) +
  geom_bar(stat="identity") + ggtitle("Cluster 0 vs. 1")
p

#FG_0vs2
log_FC_2 = c(0.650380149,
             0.579105053,
             0.351875897,
             0.302175367,
             0.294793691,
             0.277924672,
             0.270703315,
             0.237021893,
             0.220075874,
             0.185473351)

sig_name = c("C[T>G]A",
             "C[C>G]A",
             "C[C>G]C",
             "C[T>G]C",
             "A[C>G]G",
             "C[C>G]G",
             "G[C>T]G",
             "G[T>A]C",
             "A[C>A]G",
             "G[C>T]A")

FG_0v2 <- data.frame(sig_name,log_FC_2)

p <- ggplot(data=FG_0v2, aes(x=sig_name, y=log_FC_2, fill = sig_name)) +
  geom_bar(stat="identity") + ggtitle("Cluster 0 vs. 2")
p


#FG_0vs3
log_FC_2 = c(0.070561792,
             0.066665381,
             0.031712181,
             0.033222558,
             0.022273909,
             0.035439714,
             0.028657194,
             0.031836434,
             0.039898411,
             0.018685857)

sig_name = c("C[T>G]A",
             "C[C>G]A",
             "C[C>G]C",
             "C[T>G]C",
             "G[T>A]C",
             "A[C>G]G",
             "T[C>A]G",
             "T[C>G]G",
             "C[T>G]T",
             "G[C>T]A")

FG_0v3 <- data.frame(sig_name,log_FC_2)

p <- ggplot(data=FG_0v3, aes(x=sig_name, y=log_FC_2, fill = sig_name)) +
  geom_bar(stat="identity") + ggtitle("Cluster 0 vs. 3")
p

#FG_1v2
log_FC_2 = c(0.294706175,
      0.267246697,
      0.215852797,
      0.184718496,
      0.173281411,
      0.167331117,
      0.164548944,
      0.158301779,
      0.150913853,
      0.147313048)
sig_name = c("C[T>G]T",
      "C[C>G]T",
      "T[C>G]G",
      "C[T>G]A",
      "T[C>T]G",
      "C[C>A]G",
      "C[C>G]A",
      "T[C>A]G",
      "T[T>C]G",
      "C[T>C]G")
FG_1v2 <- data.frame(sig_name,log_FC_2)

p <- ggplot(data=FG_1v2, aes(x=sig_name, y=log_FC_2, fill = sig_name)) +
  geom_bar(stat="identity") + ggtitle("Cluster 1 vs. 2")
p

#FG_1vs3
log_FC_2 = c(0.162493478,
             0.120353319,
             0.111783439,
             0.095927894,
             0.088072508,
             0.087712664,
             0.085252831,
             0.078221107,
             0.069907557,
             0.067279715)

sig_name = c("C[T>G]G",
             "A[C>G]G",
             "C[C>G]G",
             "A[C>A]G",
             "G[C>G]G",
             "G[T>G]C",
             "A[T>C]G",
             "G[T>A]C",
             "A[C>T]G",
             "G[T>G]G")

FG_1v3 <- data.frame(sig_name,log_FC_2)

p <- ggplot(data=FG_1v3, aes(x=sig_name, y=log_FC_2, fill = sig_name)) +
  geom_bar(stat="identity") + ggtitle("Cluster 1 vs. 3")
p
