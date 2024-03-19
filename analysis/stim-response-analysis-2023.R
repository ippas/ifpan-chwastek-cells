setwd("/projects/ifpan-chwastek-cells")

require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)
require(ggvenn)


#save raw data:

fpkms.log <- read.table("results/log2_fpkms.csv", header = TRUE)
samples <- read.table("results/full_sample_info.csv", header = TRUE, sep = ",")

# select samples - just the OA group and just the IFNGg and LPS

samples %>% filter(
  (disease == 'O') & (treament != 'tnfa')
)-> samples

write.table(samples[,-6] , "results/selected-samples-2023.csv")

### ANALYSIS ####

results <- data.frame(fpkms.log$gene_stable_ID, fpkms.log$gene_symbol)

# 3-way anova analysis with RM:

stat <- function(counts,
                 treatment,
                 patient) {
  if (is.na(counts[1])) {
    c(NA)
  } else {
    unlist(summary(aov(as.numeric(counts) ~ 
                         as.factor(treatment)
                       + Error(
                         as.factor(patient)/as.factor(treatment))
    )
    )
    )[c(14)]
  }
  }

apply(fpkms.log[,samples$sampleID],
      1,
      stat,
      treatment = samples$treament,
      patient = samples$patient) %>% t %>% apply(.,
                                                 2,
                                                 p.adjust,
                                                 method="fdr") %>% 
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results) <- c('gene_stable_ID', 'gene_symbol', 'p.treatment')

stat.paired.t <- function(x, y) {
  if (is.na(x[1])) { NA
  } else {
    tryCatch((t.test(x ~ y))$p.value, error=function(err) NA)
  }}

ctrl_vs_lps <- which(samples$treament == "ctrl" | samples$treament == "lps")
ctrl_vs_i <- which(samples$treament == "ctrl" | samples$treament == "ifng")

fpkms.log[,samples$sampleID[ctrl_vs_lps]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=samples$treament[ctrl_vs_lps]) -> results$t.test.ctrl.vs.lps

fpkms.log[,samples$sampleID[ctrl_vs_i]] %>% data.matrix() %>%
  apply(1, stat.paired.t, y=samples$treament[ctrl_vs_i]) -> results$t.test.ctrl.vs.i


#add direction of change
check.direction <- function(x, group_1, group_2) {
  ifelse(
    mean(as.numeric(x[group_1])) -
      mean(as.numeric(x[group_2])) 
    > 0, "DOWN", "UP")
}


results$directions.ctrl.vs.lps <- apply(
  fpkms.log,
  1,
  check.direction,
  group_1 = samples$sampleID[samples$treament == 'ctrl'],
  group_2 =  samples$sampleID[samples$treament == 'lps']
)

results$directions.ctrl.vs.ifng<- apply(
  fpkms.log,
  1,
  check.direction,
  group_1 = samples$sampleID[samples$treament == 'ctrl'],
  group_2 =  samples$sampleID[samples$treament == 'ifng']
)


fold.change <- function(x, group, ctrl) {
  abs(mean(as.numeric(x[samples$sampleID[samples$treament == ctrl]]))
      -
        mean(as.numeric(x[samples$sampleID[samples$treament == group]])))
}

results %>% mutate(
  fold.change.lps=
    apply(fpkms.log, 1, fold.change, group = "lps", ctrl = "ctrl"),
  fold.change.ifng=
    apply(fpkms.log, 1, fold.change, group = "ifng", ctrl = "ctrl"),
) -> results



get.sd.max <- function(counts) {
  my.groups <- unique(samples$treament)
  sds <- vector(mode="numeric", length=3)
  for (i in seq_along(my.groups)) {
    sds[i] <- (sd(counts[samples$sampleID[which(samples$treament == my.groups[i])]]))
  }
  return(max(sds))
}

results$sd_max <- apply(fpkms.log, 1, get.sd.max)

results <- bind_cols(results, fpkms.log[,samples$sampleID])

results %>% filter(p.treatment < 0.01) -> top.genes

write.table(top.genes, './results/top-genes-2023.csv')

top.genes %>% filter(t.test.ctrl.vs.lps < 0.01) -> top.genes.lps
top.genes %>% filter(t.test.ctrl.vs.i < 0.01) -> top.genes.i

top.genes.lps %>% filter(directions.ctrl.vs.lps == 'UP') %>% select(gene_symbol) %>%
  write.table(
    .,
    file='./results/lps_up.txt',
    row.names = FALSE,
    quote = FALSE)
#%>% nrow()

top.genes.lps %>% filter(directions.ctrl.vs.lps == 'DOWN') %>% select(gene_symbol) %>%
  write.table(
    .,
    file='./results/lps_down.txt',
    row.names = FALSE,
    quote = FALSE) #%>% nrow()

top.genes.i %>% filter(directions.ctrl.vs.ifng == 'UP') %>% select(gene_symbol) %>%
  write.table(
    .,
    file='./results/ifng_up.txt',
    row.names = FALSE,
    quote = FALSE) #%>% nrow()

top.genes.i %>% filter(directions.ctrl.vs.ifng == 'DOWN') %>% select(gene_symbol) %>%
  write.table(
    .,
    file='./results/ifng_down.txt',
    row.names = FALSE,
    quote = FALSE) #%>% nrow()


### venn plotting

top.genes$lps_up <- ((top.genes$t.test.ctrl.vs.lps < 0.01) & (top.genes$directions.ctrl.vs.lps == 'UP'))
top.genes$lps_down <- ((top.genes$t.test.ctrl.vs.lps < 0.01) & (top.genes$directions.ctrl.vs.lps == 'DOWN'))

top.genes$ifng_up <- ((top.genes$t.test.ctrl.vs.i < 0.01) & (top.genes$directions.ctrl.vs.ifng == 'UP'))
top.genes$ifng_down <- ((top.genes$t.test.ctrl.vs.i < 0.01) & (top.genes$directions.ctrl.vs.ifng == 'DOWN'))

top.genes$lps <- (top.genes$t.test.ctrl.vs.lps < 0.01)
top.genes$ifng <- (top.genes$t.test.ctrl.vs.i < 0.01)

ggplot() +
  geom_venn(
    aes_string(A = 'lps_down', B =  'ifng_up'),
    data = top.genes,
    fill_color = c('#0073C2FF', '#EFC000FF'),
   set_names = c('regulated by LPS', 'regulated by IFNg'),
    #set_name_size = 4,
    auto_scale = TRUE,
    show_percentage = FALSE,
    show_outside = 'none',
  ) + theme_void()


### HEATMAP PLOTTING 

top.genes %>% filter((t.test.ctrl.vs.lps < 0.00000001) | (t.test.ctrl.vs.i < 0.00000001)) %>%
  filter(
    (
      (fold.change.lps > 3.5) | (fold.change.lps < -3.5)
      ) | (
        (fold.change.ifng > 3.5) | (fold.change.ifng < -3.5)
      )
  ) -> selected
  

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- c("CTRL", "IFNg", "LPS")

col.labels <- c(rep("", 10), group.names[1], rep(" ", 14), 
                group.names[2], rep(" ", 14),
                group.names[3], rep(" ", 10)
                )


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}

selected[,samples$sampleID[order(samples$treament)]] %>% data.matrix() -> to.plot
rownames(to.plot) <- coalesce(na_if(selected$gene_symbol, ""), selected$gene_stable_ID)

to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 2.5) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(16,32),
    sepwidth = c(0.3,0.3),
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.6,
    offsetCol = 0,
    cexCol = 1
  )


### read in a gene list of secreted factors
secreted <- read.table('./data/secreted_genes.csv')

secreted$top <- secreted$V1 %in% top.genes$gene_symbol
secreted$secreted <- TRUE

ggplot() +
  geom_venn(
    aes_string(A = 'secreted', B =  'top'),
    data = secreted,
    fill_color = c('#0073C2FF', '#EFC000FF'),
    set_names = c('genes with secreted protein protucts', 'regulated by LPS and/or IFNg'),
    #set_name_size = 4,
    auto_scale = TRUE,
    show_percentage = FALSE,
    show_outside = 'none',
  ) + theme_void()

fpkms.log %>% select(gene_symbol) %>% unique()
