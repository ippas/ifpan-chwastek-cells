setwd("/projects/ifpan-chwastek-cells")

require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)

# create sample info:
samples_cuffnorm <- data.frame(read.table('data/cuffnorm/samples.table', header = TRUE))

samples_cuffnorm$file %>% as.character() %>%
  strsplit("[/]") %>% data.frame() %>% t() -> samples_temp

samples_temp[,4] %>% strsplit("[.]") %>% data.frame() %>% t() -> samples_temp
samples_cuffnorm$sample <- samples_temp[,1]
rm(samples_temp)

samples <- data.frame(samples_cuffnorm$sample_id, samples_cuffnorm$sample)
rm(samples_cuffnorm)

#assign groups

colnames(samples) <- c('cuffnormID', 'sampleID')

samples$disease <- substr(samples$sampleID,1,1)
samples$patient <- substr(samples$sampleID,2,str_length(samples$sampleID)-1)
samples$treament <- substr(samples$sampleID,str_length(samples$sampleID),str_length(samples$sampleID))

samples$treament[samples$treament == 1 ] <- 'ctrl'
samples$treament[samples$treament == 2 ] <- 'ifng'
samples$treament[samples$treament == 3 ] <- 'tnfa'
samples$treament[samples$treament == 4 ] <- 'lps'

# switch wrongly named samples: O662 i O664, O1212 i O1214, O1222 i O1224, O1842 i O1844

samples$wrong_sampleID <- samples$sampleID

switch_sample_ids <- function(x,y,z) {
  z <- as.character(z)
  z[z == x] <- "A"
  z[z == y] <- "B"
  z[z == "A"] <- y
  z[z == "B"] <- x
  return(z)
}

samples$sampleID <- switch_sample_ids("O662", "O664", samples$sampleID) 
samples$sampleID <- switch_sample_ids("O1212", "O1214", samples$sampleID)
samples$sampleID <- switch_sample_ids("O1222", "O1224", samples$sampleID)
samples$sampleID <- switch_sample_ids("O1842", "O1844", samples$sampleID)


# read reads 
genes <- data.frame(read.table('data/biomart-human.txt', sep = "\t", header = TRUE))

fpkms <- data.frame(read.table('data/cuffnorm/genes.fpkm_table', header = TRUE))

fpkms <- cbind(fpkms$tracking_id, genes$HGNC.symbol[match(fpkms$tracking_id, genes$Gene.stable.ID)],fpkms[,as.character(samples$cuffnormID)])
colnames(fpkms)[c(1,2)] <- c('gene_stable_ID', 'gene_symbol')
colnames(fpkms)[c(3:106)] <- as.character(samples$sampleID[match(colnames(fpkms)[c(3:106)], samples$cuffnormID)])

fpkms.normalised <- data.matrix(fpkms[,c(-1,-2)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA

fpkms.log <- data.frame(as.character(fpkms$gene_stable_ID), as.character(fpkms$gene_symbol), data.matrix(fpkms.log))

colnames(fpkms.log) <- colnames(fpkms)


#save raw data:
write.table(fpkms, "results/raw_fpkms.csv", row.names = FALSE)
write.table(fpkms.log, "results/log2_fpkms.csv", row.names = FALSE)
write.csv(samples, "results/full_sample_info.csv", row.names = FALSE)


### ANALYSIS ####

results <- data.frame(fpkms.log$gene_stable_ID, fpkms.log$gene_symbol)

# 3-way anova analysis with RM:

stat <- function(counts,
                 treatment,
                 disease,
                 patient) {
  if (is.na(counts[1])) {
    c(rep(NA, 3))
  } else {
    unlist(summary(aov(counts ~ 
                         as.factor(treatment)
                       * as.factor(disease)
                       + Error(as.factor(patient)/as.factor(treatment)))))[c(9,23,24)]}}

apply(fpkms.log[,as.character(samples$sampleID)],
      1,
      stat,
      treatment = samples$treament,
      disease = samples$disease,
      patient = samples$patient) %>% t %>% apply(.,
                                                 2,
                                                 p.adjust,
                                                 method="fdr") %>% 
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results) <- c('gene_stable_ID', 'gene_symbol', 'p.disease','p.treatment','p.interaction')

# question 1: difference between patients and controls
samples$group_full <- paste(samples$treament, "_", samples$disease, sep="")

chosen_samples <- as.factor(samples$group_full[which(samples$group_full == "ctrl_H" | samples$group_full == "ctrl_O")])


stat.paired.t <- function(x, y) {
  if (is.na(x[1])) { NA
  } else {
    tryCatch((t.test(x ~ y))$p.value, error=function(err) NA)
  }}


fpkms.log[,samples$sampleID[which(samples$group_full == "ctrl_H" | samples$group_full == "ctrl_O")]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=chosen_samples) -> results$t.test.ctrl.h.vs.o

#add direction change beween controls

check.direction <- function(x) {
  ifelse(
    mean(as.numeric(x[samples$sampleID[samples$group_full == "ctrl_H"]])) -
      mean(as.numeric(x[samples$sampleID[samples$group_full == "ctrl_O"]])) 
    > 0, "DOWN", "UP")
}


results$directions.ctrl.h.vs.o <- apply(fpkms.log, 1, check.direction)


# question 2: different responses between treatments:
# for this analysis we will look at difference to fold changes for different treatments
# calculate fold changes:

fold.change.special <- function(x, group, ctrl) {
  abs(mean(as.numeric(x[samples$sampleID[samples$group_full == ctrl]]))
      -
        mean(as.numeric(x[samples$sampleID[samples$group_full == group]])))
}

results %>% mutate(
  fold.change.tnf.h=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_H", ctrl = "ctrl_H"),
  fold.change.tnf.o=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_O", ctrl = "ctrl_O"),
  fold.change.lps.h=
    apply(fpkms.log, 1, fold.change.special, group = "lps_H", ctrl = "ctrl_H"),
  fold.change.lps.o=
    apply(fpkms.log, 1, fold.change.special, group = "lps_O", ctrl = "ctrl_O"),
  fold.change.ifg.h=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_H", ctrl = "ctrl_H"),
  fold.change.ifg.o=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_O", ctrl = "ctrl_O"),
  fold.change.tnf.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_H", ctrl = "tnfa_O"),
  fold.change.lps.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "lps_H", ctrl = "lps_O"),
  fold.change.ifg.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_H", ctrl = "ifng_O")
) -> results


#calculate differences in fold changes for each treatment:

results %>% mutate(
  difference.tnf= abs(fold.change.tnf.h - fold.change.tnf.o),
  difference.ifg= abs(fold.change.ifg.h - fold.change.ifg.o),
  difference.lps= abs(fold.change.lps.h - fold.change.lps.o)
  ) -> results


selected.genes <- c('CXCL16', 'CCL4', 'IGF2', 'CXCL6', 'CXCL10', 'TNFSF10', 'SRGN', 'DKK1')


get.sd.max <- function(counts) {
  my.groups <- unique(samples$group_full)
  sds <- vector(mode="numeric", length=8)
  for (i in seq_along(my.groups)) {
    sds[i] <- (sd(counts[samples$sampleID[which(samples$group_full == my.groups[i])]]))
  }
  return(max(sds))
}

results$sd <- apply(fpkms.log, 1, get.sd.max)

results <- bind_cols(results, fpkms.log[,c(3:106)])

results %>% filter(gene_symbol %in% selected.genes) -> selected.results
results %>% filter((p.disease < 0.3) & (p.treatment < 0.05) & (fold.change.tnf.h.vs.o > 1) & (p.interaction < 0.9) & (sd < 2.75)) -> top.genes


write.table(results, "results/oa_cells_results_full.csv", row.names = FALSE)
write.table(selected.results, "results/selected_oa_results_full.csv", row.names = FALSE)
write.table(top.genes, "results/selected_oa_150_results_full.csv", row.names = FALSE)


### at this point IFNg results are removed from the analysis ###

to.drop <- c("fold.change.ifg.h", 
             "fold.change.ifg.o", 
             "fold.change.ifg.h.vs.o",
             "difference.ifg",
             samples$sampleID[which((samples$group_full == 'ifng_H') | (samples$group_full == 'ifng_O'))])

results <- results[ , !(names(results) %in% to.drop)]

results %>% filter(gene_symbol %in% selected.genes) -> selected.results
results %>% filter((p.disease < 0.3) & (p.treatment < 0.05) & (fold.change.tnf.h.vs.o > 1) & (p.interaction < 0.9) & (sd < 2.75)) -> top.genes

write.table(results, "results/oa_cells_results.csv", row.names = FALSE)
write.table(selected.results, "results/selected_oa_results.csv", row.names = FALSE)
write.table(top.genes, "results/selected_oa_150_results.csv", row.names = FALSE)


### HEATMAP PLOTTING (this is done without ifng)

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- c("ctrl_H", "ctrl_O", "lps_H",  "lps_O",  "tnfa_H", "tnfa_O")

col.labels <- c(rep("", 5), group.names[1], rep(" ", 12), 
                group.names[2], rep(" ", 14),
                group.names[3], rep(" ", 10),
                group.names[4], rep(" ", 12),
                group.names[5], rep(" ", 12),
                group.names[6], rep(" ", 12))


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}

samples.short <- samples[(samples$group_full != 'ifng_H') & (samples$group_full != 'ifng_O'), ]

fpkms.log[match(
  results$gene_stable_ID[which((results$p.disease < 0.3) & (results$p.treatment < 0.05) & (results$fold.change.tnf.h.vs.o > 1) & (results$p.interaction <0.9) & (results$sd < 2.75))],
  fpkms.log$gene_stable_ID),
  as.character(samples.short$sampleID[order(as.character(samples.short$group_full))])] %>% data.matrix()-> to.plot


to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(10,26,36,52,62,78,88),
    sepwidth = c(0.3,0.3),
    labRow=fpkms.log$gene_symbol[match(rownames(.), rownames(fpkms.log))],
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.3,
    offsetCol = 0
  )


#### rerun the analyses with smaller groups ####

#OA (22, 23, 66, 88, 117) H (09, 08, 05, 06, 07)

samples.n5 <- samples[samples$patient %in% c("09", "08", "05", "06", "07", "22", "23",
                                             "66", "88", "117"),]

results.n5 <- data.frame(fpkms.log$gene_stable_ID, fpkms.log$gene_symbol)

apply(fpkms.log[,as.character(samples.n5$sampleID)],
      1,
      stat,
      treatment = samples.n5$treament,
      disease = samples.n5$disease,
      patient = samples.n5$patient) %>% t %>% apply(.,
                                                 2,
                                                 p.adjust,
                                                 method="fdr") %>% 
  data.frame() %>%
  bind_cols(results.n5,.) -> results.n5


colnames(results.n5) <- c('gene_stable_ID', 'gene_symbol', 'p.disease','p.treatment','p.interaction')

# question 1: difference between patients and controls

chosen_samples <- as.factor(samples.n5$group_full[which(samples.n5$group_full == "ctrl_H" | samples.n5$group_full == "ctrl_O")])


fpkms.log[,samples.n5$sampleID[which(samples.n5$group_full == "ctrl_H" | samples.n5$group_full == "ctrl_O")]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=chosen_samples) -> results.n5$t.test.ctrl.h.vs.o

#add direction change beween controls

check.direction <- function(x) {
  ifelse(
    mean(as.numeric(x[samples.n5$sampleID[samples.n5$group_full == "ctrl_H"]])) -
      mean(as.numeric(x[samples.n5$sampleID[samples.n5$group_full == "ctrl_O"]])) 
    > 0, "DOWN", "UP")
}


results.n5$directions.ctrl.h.vs.o <- apply(fpkms.log, 1, check.direction)


# question 2: different responses between treatments:
# for this analysis we will look at difference to fold changes for different treatments
# calculate fold changes:

fold.change.special <- function(x, group, ctrl) {
  abs(mean(as.numeric(x[samples.n5$sampleID[samples.n5$group_full == ctrl]]))
      -
        mean(as.numeric(x[samples.n5$sampleID[samples.n5$group_full == group]])))
}

results.n5 %>% mutate(
  fold.change.tnf.h=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_H", ctrl = "ctrl_H"),
  fold.change.tnf.o=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_O", ctrl = "ctrl_O"),
  fold.change.lps.h=
    apply(fpkms.log, 1, fold.change.special, group = "lps_H", ctrl = "ctrl_H"),
  fold.change.lps.o=
    apply(fpkms.log, 1, fold.change.special, group = "lps_O", ctrl = "ctrl_O"),
  fold.change.ifg.h=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_H", ctrl = "ctrl_H"),
  fold.change.ifg.o=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_O", ctrl = "ctrl_O"),
  fold.change.tnf.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "tnfa_H", ctrl = "tnfa_O"),
  fold.change.lps.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "lps_H", ctrl = "lps_O"),
  fold.change.ifg.h.vs.o=
    apply(fpkms.log, 1, fold.change.special, group = "ifng_H", ctrl = "ifng_O")
) -> results.n5


#calculate differences in fold changes for each treatment:

results.n5 %>% mutate(
  difference.tnf= abs(fold.change.tnf.h - fold.change.tnf.o),
  difference.ifg= abs(fold.change.ifg.h - fold.change.ifg.o),
  difference.lps= abs(fold.change.lps.h - fold.change.lps.o)
) -> results.n5


get.sd.max <- function(counts) {
  my.groups <- unique(samples.n5$group_full)
  sds <- vector(mode="numeric", length=8)
  for (i in seq_along(my.groups)) {
    sds[i] <- (sd(counts[samples.n5$sampleID[which(samples.n5$group_full == my.groups[i])]]))
  }
  return(max(sds))
}

results.n5$sd <- apply(fpkms.log, 1, get.sd.max)

results.n5 <- bind_cols(results.n5, fpkms.log[,samples.n5$sampleID])

results.n5 %>% filter(gene_symbol %in% selected.genes) -> selected.results.n5
results.n5 %>% filter((fold.change.tnf.h.vs.o > 0.95) & (p.treatment < 0.4 ) & (p.disease < 0.75) & (sd < 3.5)) -> top.genes.n5

write.table(results.n5, "results/oa_cells_results_full_n5.csv", row.names = FALSE)
write.table(selected.results.n5, "results/selected_oa_results_full_n5.csv", row.names = FALSE)
write.table(top.genes.n5, "results/selected_oa_150_results_full_n5.csv", row.names = FALSE)


### at this point IFNg results are removed from the analysis ###

to.drop <- c("fold.change.ifg.h", 
             "fold.change.ifg.o", 
             "fold.change.ifg.h.vs.o",
             "difference.ifg",
             samples.n5$sampleID[which((samples.n5$group_full == 'ifng_H') | (samples.n5$group_full == 'ifng_O'))])

results.n5 <- results.n5[ ,!(names(results.n5) %in% to.drop)]

results.n5 %>% filter(gene_symbol %in% selected.genes) -> selected.results.n5
results.n5 %>% filter((fold.change.tnf.h.vs.o > 0.95) & (p.treatment < 0.4 ) & (p.disease < 0.75) & (sd < 3.5)) -> top.genes.n5

write.table(results.n5, "results/oa_cells_results_n5.csv", row.names = FALSE)
write.table(selected.results.n5, "results/selected_oa_results_n5.csv", row.names = FALSE)
write.table(top.genes.n5, "results/selected_oa_150_results_n5.csv", row.names = FALSE)

### HEATMAP PLOTTING for n5(this is done without ifng)

samples.short.n5 <- samples.n5[(samples.n5$group_full != 'ifng_H') & (samples.n5$group_full != 'ifng_O'), ]

fpkms.log[match(
  results.n5$gene_stable_ID[which((results.n5$p.disease < 0.75) & (results.n5$p.treatment < 0.4) & (results.n5$fold.change.tnf.h.vs.o > 0.95) & (results$sd < 3.5))],
  fpkms.log$gene_stable_ID),
  as.character(samples.short.n5$sampleID[order(as.character(samples.short.n5$group_full))])] %>% data.matrix()-> to.plot

group.names <- c("ctrl_H", "ctrl_O", "lps_H",  "lps_O",  "tnfa_H", "tnfa_O")

col.labels <- c(rep("", 2), group.names[1], rep(" ", 4), 
                group.names[2], rep(" ", 4),
                group.names[3], rep(" ", 4),
                group.names[4], rep(" ", 4),
                group.names[5], rep(" ", 4),
                group.names[6], rep(" ", 3))


to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(5,10,15,20,25,30,35),
    sepwidth = c(0.3,0.3),
    labRow=fpkms.log$gene_symbol[match(rownames(.), rownames(fpkms.log))],
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.2,
    offsetCol = 0
  )


