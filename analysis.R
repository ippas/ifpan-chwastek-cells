
require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)

# create sample info:
samples_cuffnorm <- data.frame(read.table('samples.table', header = TRUE))
samples <- data.frame(read.table('samples.csv', header = TRUE))

samples_cuffnorm$file %>% as.character() %>%
  strsplit("[/]") %>% data.frame() %>% t() -> samples_temp

samples_temp[,4] %>% strsplit("[.]") %>% data.frame() %>% t() -> samples_temp
samples_cuffnorm$sample <- samples_temp[,1]
rm(samples_temp)

samples$cuffnorm_id <- samples_cuffnorm$sample_id[match(samples$SampleID, samples_cuffnorm$sample)]
rm(samples_cuffnorm)

colnames(samples) <- c('sampleID', 'group', 'cuffnormID')

samples$disease <- substr(samples$sampleID,1,1)
samples$patient <- substr(samples$sampleID,2,str_length(samples$sampleID)-1)



genes <- data.frame(read.table('biomart.txt', sep = "\t", header = TRUE))

fpkms <- data.frame(read.table('genes.fpkm_table', header = TRUE))

fpkms <- cbind(fpkms$tracking_id, genes$HGNC.symbol[match(fpkms$tracking_id, genes$Gene.stable.ID)],fpkms[,as.character(samples$cuffnormID)])

colnames(fpkms)[c(1,2)] <- c('gene_stable_ID', 'gene_symbol')
colnames(fpkms)[c(3:106)] <- as.character(samples$sampleID[match(colnames(fpkms)[c(3:106)], samples$cuffnormID)])

fpkms.normalised <- data.matrix(fpkms[,c(-1,-2)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised,copy=FALSE)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA

fpkms.log <- data.frame(cbind(as.character(fpkms$gene_stable_ID), as.character(fpkms$gene_symbol), fpkms.log))

colnames(fpkms.log)[c(1,2)] <- c('gene_stable_ID', 'gene_symbol')


# replace wrongly named  in sample info samples: O662 i O664, O1212 i O1214, O1222 i O1224, O1842 i O1844

samples$previous_sampleID <- samples$sampleID

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

#save raw data:
write.table(fpkms, "raw_fpkms.csv", row.names = FALSE)
write.table(fpkms.log, "log2_fpkms.csv", row.names = FALSE)
write.csv(samples, "full_sample_info.csv", row.names = FALSE)


#exploratory analysis:

samples$group_full <- paste(samples$group, "_", samples$disease, sep="")

stat.one.way <- function(counts, groups) {
  ifelse(
    is.na(counts[1]),
    NA,
    unlist(summary(aov(counts ~ as.factor(groups))))[9]) }


results <- data.frame(fpkms.log$gene_stable_ID, fpkms.log$gene_symbol)

#perform one-way stat on all groups:

apply(fpkms.log[,match(samples$sampleID, colnames(fpkms.log))],
      1,
      stat.one.way,
      groups=samples$group_full) %>% 
  p.adjust(method="fdr") -> results$p.one.way



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
      treatment = samples$group,
      disease = samples$disease,
      patient = samples$patient) %>% t %>% apply(.,
                  2,
                  p.adjust,
                  method="fdr") %>% 
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results)[c(4:6)] <- c('p.disease','p.treatment','p.interaction')

# question 1: difference between patients and controls

chosen_samples <- as.factor(samples$group_full[which(samples$group_full == "Ctrl_H" | samples$group_full == "Ctrl_O")])


stat.paired.t <- function(x, y) {
  if (is.na(x[1])) { NA
  } else {
    tryCatch((t.test(x ~ y))$p.value, error=function(err) NA)
  }}


fpkms.log[,samples$sampleID[which(samples$group_full == "Ctrl_H" | samples$group_full == "Ctrl_O")]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=chosen_samples) -> results$t.test.ctrl.h.vs.o

#add direction change beween controls

check.direction <- function(x) {
ifelse(
  mean(as.numeric(x[samples$sampleID[samples$group_full == "Ctrl_H"]])) -
  mean(as.numeric(x[samples$sampleID[samples$group_full == "Ctrl_O"]])) 
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
    apply(fpkms.log, 1, fold.change.special, group = "TNFa_H", ctrl = "Ctrl_H"),
  fold.change.tnf.o=
    apply(fpkms.log, 1, fold.change.special, group = "TNFa_O", ctrl = "Ctrl_O"),
  fold.change.lps.h=
    apply(fpkms.log, 1, fold.change.special, group = "LPS_H", ctrl = "Ctrl_H"),
  fold.change.lps.o=
    apply(fpkms.log, 1, fold.change.special, group = "LPS_O", ctrl = "Ctrl_O"),
  fold.change.ifg.h=
    apply(fpkms.log, 1, fold.change.special, group = "IFg_H", ctrl = "Ctrl_H"),
  fold.change.ifg.o=
    apply(fpkms.log, 1, fold.change.special, group = "IFg_O", ctrl = "Ctrl_O")
) -> results


#calculate differences in fold changes for each treatment:

results %>% mutate(
  difference.tnf= abs(fold.change.tnf.h - fold.change.tnf.o),
  difference.ifg= abs(fold.change.ifg.h - fold.change.ifg.o),
  difference.lps= abs(fold.change.lps.h - fold.change.lps.o)
  ) -> results


results <- bind_cols(results, fpkms.log[,c(3:106)])
write.table(results, "oa_cells_results.csv", row.names = FALSE)

### select top gene list:

results %>% 
  filter(p.disease < 0.1 & t.test.ctrl.h.vs.o < 0.0001 & directions.ctrl.h.vs.o == "UP") %>% 
  select(fpkms.log.gene_symbol) %>% write.table(row.names = FALSE, quote = FALSE)

results %>% 
  filter(p.disease < 0.1 & t.test.ctrl.h.vs.o < 0.0001 & directions.ctrl.h.vs.o == "DOWN") %>% 
  select(fpkms.log.gene_symbol) %>% write.table(row.names = FALSE, quote = FALSE)

results %>% 
  filter(p.disease < 0.1 & t.test.ctrl.h.vs.o < 0.0001) %>% 
  write.table("healthy_vs_disease.csv", row.names = FALSE)


results %>% 
  filter(p.interaction < 0.1) %>% 
  write.table("top_interaction.csv", row.names = FALSE)



### HEATMAP PLOTTING:

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(samples$group_full[order(samples$group_full)])



col.labels <- c(rep("", 5), group.names[1], rep(" ", 12), 
                group.names[2], rep(" ", 14),
                group.names[3], rep(" ", 10),
                group.names[4], rep(" ", 12),
                group.names[5], rep(" ", 12),
                group.names[6], rep(" ", 12),
                group.names[7], rep(" ", 12),
                group.names[8], rep(" ", 12))

for (i in unique(samples$group_full[order(samples$group_full)])) {
  print(i)
  print(sum(str_count(samples$group_full, i)))
}



cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.one.way < 0.00000000000000003)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot


fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.interaction < 0.1 & results$difference.ifg > 0.75)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot


fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.interaction < 0.1 & results$difference.tnf > 0.75)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot

fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.interaction < 0.1 & results$difference.lps > 0.75)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot

fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.interaction < 0.1)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot




fpkms.log[match(
  results$fpkms.log.gene_stable_ID[which(results$p.disease < 0.1 & results$t.test.ctrl.h.vs.o < 0.0001)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot


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
    offsetCol = 0.1
  )



