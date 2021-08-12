### this is an outdated filed that depends on data no longer available for archive purposes ###

require(magrittr)
require(dplyr)
require(lmerTest)
require(RColorBrewer)
require(gplots)

# a heatmap of a list of genes from Kuba

results <- data.frame(read.table("oa_cells_results.csv", header = TRUE))
fpkms.log <- data.frame(read.table("log2_fpkms.csv", header = TRUE))

genes.1 <- data.frame(read.table('genes-1.csv', header = FALSE))
genes.2 <- data.frame(read.table('genes-2.csv', header = FALSE))
  
results %>% filter(fpkms.log$gene_symbol %in% genes.1$V1) -> selected.results.1
results %>% filter(fpkms.log$gene_symbol %in% genes.2$V1) -> selected.results.2

#load sample data
samples <- data.frame(read.table("full_sample_info.csv", header = TRUE, sep = ","))
pheno <- data.frame(read.table('patient_data.csv', header = TRUE))
pheno$sample <- gsub("[^0-9]", "", pheno$sample)
samples$sex <-pheno$gender[match(as.numeric(samples$patient), as.numeric(pheno$sample))]
samples$age <-as.numeric(pheno$gender.1)[match(as.numeric(samples$patient), as.numeric(pheno$sample))]
samples$hemoglobin <-pheno$Hemoglobin[match(as.numeric(samples$patient), as.numeric(pheno$sample))]
samples$leukocytes <-pheno$Leukocyte[match(as.numeric(samples$patient), as.numeric(pheno$sample))]
samples$group_full <- paste(samples$group, "_", samples$disease, sep="")

# run stats
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


apply(selected.results.1[,as.character(samples$sampleID)],
      1,
      full_model) %>% t %>% apply(.,
                                  2,
                                  p.adjust,
                                  method="fdr") %>% 
  data.frame() -> temp.full.model.selected


apply(selected.results.2[,as.character(samples$sampleID)],
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


#cbind(selected.results.1, temp.full.model.selected) -> selected.results.1

cbind(selected.results.2, temp.full.model.selected) -> selected.results.2


fpkms.log[match(
  selected.results.1$fpkms.log.gene_stable_ID[which(selected.results.1$full.model.fdr.disease.treatment < 0.05)],
  fpkms.log$gene_stable_ID),as.character(samples$sampleID[order(as.character(samples$group_full))])] %>%
  na.omit() %>% data.matrix()-> to.plot


fpkms.log[match(
  selected.results.2$fpkms.log.gene_stable_ID[which(selected.results.2$full.model.fdr.disease.treatment < 0.05)],
  fpkms.log$gene_stable_ID),as.character(samples$sampleID[order(as.character(samples$group_full))])] %>%
  na.omit() %>% data.matrix()-> to.plot



#plot_results:

### HEATMAP PLOTTING:
mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)


greys_scale <- brewer.pal(9,"Greys")
more_greys <- colorRampPalette(greys_scale)

group.names <- unique(samples$group_full[order(samples$group_full)])


col.labels <- c(rep("", 5), group.names[1], rep(" ", 12), 
                group.names[2], rep(" ", 14),
                group.names[3], rep(" ", 10),
                group.names[4], rep(" ", 12),
                group.names[5], rep(" ", 12),
                group.names[6], rep(" ", 12),
                group.names[7], rep(" ", 12),
                group.names[8], rep(" ", 12))




cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


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
    cexRow = 0.8,
    cexCol = 0.8,
    offsetCol = 0.1
  )



full.model.all <- data.frame(read.table("full_model_export.csv", header = TRUE))
