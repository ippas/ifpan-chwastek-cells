### this is an outdated filed that depends on some data that is no longer available and is kept for archive purposes ###


require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(dplyr)
require(stringi)
require(stringr)
require(lmerTest)
require(enrichR)


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

results <- data.frame(fpkms.log$gene_stable_ID, fpkms.log$gene_symbol)



# question 1: difference between patients and controls - simple t-test

chosen_samples <- as.factor(samples$group_full[which(samples$group_full == "Ctrl_H" | samples$group_full == "Ctrl_O")])


stat.paired.t <- function(x, y) {
  if (is.na(x[1])) { NA
  } else {
    tryCatch((t.test(x ~ y))$p.value, error=function(err) NA)
  }}


fpkms.log[,samples$sampleID[which(samples$group_full == "Ctrl_H" | samples$group_full == "Ctrl_O")]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=chosen_samples) -> results$t.test.ctrl.h.vs.o

#add direction change between controls

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



#add more phenotypic data

pheno <- data.frame(read.table('patient_data.csv', header = TRUE))

pheno$sample <- gsub("[^0-9]", "", pheno$sample)

samples$sex <-as.numeric(pheno$gender)[match(samples$patient, pheno$sample)]
samples$age <-as.numeric(pheno$gender.1)[match(samples$patient, pheno$sample)]
samples$hemoglobin <-pheno$Hemoglobin[match(samples$patient, pheno$sample)]
samples$leukocytes <-pheno$Leukocyte[match(samples$patient, pheno$sample)]



correlate_with <- function(x, y) { 
  ifelse(
    is.na(x[1]),
    NA,
    (unlist(cor.test(as.numeric(x), as.numeric(y))))["p.value"])
     }



apply(fpkms.log[,as.character(samples$sampleID)],
      1,
      correlate_with,
      y = samples$age) %>% 
      t %>% apply(., 2, p.adjust, method="fdr")  -> results$age.correlation

apply(fpkms.log[,as.character(samples$sampleID)],
      1,
      correlate_with,
      y = samples$sex) %>% 
  t %>% apply(., 2, p.adjust, method="fdr")  -> results$sex.correlation



### fit a GLM including age and sex as factors


full_model <- function(counts) {
  if (is.na(counts[1])) {
    c(rep(NA, 7))
  } else {
    data <- bind_cols(counts = as.numeric(counts),
                      treatment = as.factor(samples$group),
                      disease = as.factor(samples$disease),
                      sex = samples$sex,
                      age = samples$age,
                      patient = as.factor(samples$patient))
    unlist(anova(lmer(counts ~ treatment 
               + disease
               + sex
               + age
               + disease:sex
               + disease:age
               + disease:treatment
               + (1|patient), data = data)))[c("Pr(>F)1",
                                               "Pr(>F)2",
                                               "Pr(>F)3",
                                               "Pr(>F)4",
                                               "Pr(>F)5","Pr(>F)6","Pr(>F)7")]}}

full_model(fpkms.log[1,as.character(samples$sampleID)])
                
                           
apply(fpkms.log[,as.character(samples$sampleID)],
      1,
      full_model) %>% t %>% apply(.,
                                                 2,
                                                 p.adjust,
                                                 method="fdr") %>% 
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results)

colnames(results)[c(16:22)]
 


colnames(results)[c(16:22)] <- c("full.model.fdr.treatment", "full.model.fdr.disease", "full.model.fdr.sex",
                                 "full.model.fdr.age", 
                                 "full.model.fdr.disease.sex", "full.model.fdr.disease.age", "full.model.fdr.disease.treatment")

results %>% filter(full.model.fdr.disease.treatment < 0.1 & difference.tnf > 0.5) %>% nrow
results %>% filter(full.model.fdr.disease.treatment < 0.1 & difference.ifg > 0.8) %>% nrow
results %>% filter(full.model.fdr.disease.treatment < 0.1 & difference.lps > 0.7) %>% nrow

colnames(results)
#export these results just for genes with significant interaction:

results %>% filter(full.model.fdr.treatment.disease < 0.1) -> export
export <- export[,-c(4,5,6)]
export <- bind_cols(export, fpkms.log[match(export$fpkms.log.gene_stable_ID, fpkms.log$gene_stable_ID),c(3:106)])

write.table(export, "full_model_export.csv", row.names = FALSE)

### HEATMAP PLOTTING:

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)


greys_scale <- brewer.pal(9,"Greys")
more_greys <- colorRampPalette(greys_scale)


samples <- data.frame(read.table("full_sample_info.csv", header = TRUE, sep = ","))

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
  results$fpkms.log.gene_stable_ID[which(results$full.model.fdr.disease.treatment < 0.05)],
  fpkms.log$gene_stable_ID),
  as.character(samples$sampleID[order(as.character(samples$group_full))])] %>% data.matrix()-> to.plot

#plot heatmap for each of the factors:
top.results <- data.frame(read.table("top_genes.csv", header = TRUE))
fpkms.log <- data.frame(read.table("log2_fpkms.csv", header = TRUE))


#add extra filter to the results
x_larger_than_n <- function(x,n) {
(x > n) %>% as.numeric() %>% sum() -> a                                                                                                                                                            
return(a)}

top.results[,as.character(samples$sampleID)] %>% apply(1, x_larger_than_n, n=6) -> top.results$abundance.filter

write.table(top.results, "top_genes.csv", row.names = FALSE)


fpkms.log[match(
  top.results$fpkms.log.gene_stable_ID[which(top.results$difference.tnf > 0.75)],
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
    sepwidth = c(0.2,0.2),
    labRow=fpkms.log$gene_symbol[match(rownames(.), rownames(fpkms.log))],
    labCol=col.labels,
    # ColSideColors = more_greys(30)[as.factor(samples$age[order(as.character(samples$group_full))])],
    srtCol = 45,
    cexRow = 0.4,
    cexCol = 0.8,
    offsetCol = 0.1
  )



# treatment, age_group, treatment x age group


results <- bind_cols(results, fpkms.log[,c(3:106)])

write.table(results, "oa_cells_results.csv", row.names = FALSE)

results %>% filter(full.model.fdr.disease.treatment < 0.1) -> top_genes
write.table(top_genes, "top_genes.csv", row.names = FALSE)



write.table(colnames(results), row.names = FALSE, quote=FALSE)



# a heatmap of a list of genes from Kuba

results <- data.frame(read.table("oa_cells_results.csv", header = TRUE))

genes.from.kuba <- data.frame(read.table('chwastek-chosen-genes.csv', header = FALSE))

results %>% filter(fpkms.log.gene_symbol %in%  genes.from.kuba$V1) -> selected.results

full_model(selected.results[1,as.character(samples$sampleID)])


apply(selected.results[,as.character(samples$sampleID)],
      1,
      full_model) %>% t %>% apply(.,
                                  2,
                                  p.adjust,
                                  method="fdr") %>% 
  data.frame() -> temp.full.model.selected



colnames(temp.full.model.selected)

colnames(temp.full.model.selected) <- c("full.model.fdr.treatment", "full.model.fdr.disease", "full.model.fdr.sex",
                                 "full.model.fdr.age", 
                                 "full.model.fdr.disease.sex", "full.model.fdr.disease.age", "full.model.fdr.disease.treatment")


selected.results[,c(18:24)] <- temp.full.model.selected


fpkms.log[match(
  selected.results$fpkms.log.gene_stable_ID[which(selected.results$full.model.fdr.disease.treatment < 0.1)],
  fpkms.log$gene_stable_ID),as.character(samples$sampleID[order(as.character(samples$group_full))])] %>%
na.omit() %>% data.matrix()-> to.plot


# get clusters of top genes based on correlation:

top.results <- data.frame(read.table("top_genes.csv", header = TRUE))

top.results[,as.character(samples$sampleID[order(as.character(samples$group_full))])] %>%
  na.omit() %>% data.matrix()-> to.plot


#get 6 clusters:

to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  {1 - cor(.)} %>% 
  t %>%
  as.dist %>% 
  hclust %>%
  cutree(., k=6) %>% 
  as.table() %>% 
  as.data.frame() -> my.clusters

my.clusters$gene.symbol <- top.results$fpkms.log.gene_symbol

write.table(my.clusters[,c(-1)], "6_clusters_top_genes.csv", row.names = FALSE)

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
    #labRow=my.clusters$gene.symbol,
    labRow=my.clusters$Freq,
    labCol=col.labels,
    RowSideColors = more_greys(8)[c(2:8)][my.clusters$Freq],
    srtCol = 45,
    cexRow = 0.7,
    offsetCol = 0.1)


write.table(my.clusters$gene.symbol[my.clusters$Freq == 6], row.names = FALSE, quote = FALSE)
