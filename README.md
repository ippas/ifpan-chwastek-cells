## RNA-Seq in Ostheoarthritis model

Patients cells treated with 3 compounds (4 h, 10 ng/ml), 2-TNFa, 3-LPS, 4-IFg.

### Methods:

All samples were checked for quality with fastQC v0.11.9 and aligned to a GRCh38 *Homo sapiens* reference genome (release 100 from Ensembl database) with hisat2 2.2.0. Cufflinks v 2.2.1 package and GTF from the Ensembl gene database were used to quantify (cuffquant) and normalize (cuffnorm) transcripts to fpkms (Fragments Per Kilobase of transcript per Million fragments mapped).

All statisical analyses were performed with R software v3.6. [Click here for the analysis code](analysis.R). R was run in this docker container: `rocker/rstudio:3.6.3-ubuntu18.04` 

See [this file](pipeline_run_and_inputs.md) for details on preprocessing.


*note: a few samples were mislabeled in this experiment. The labels were corrected by switching sample names in sample info. Full sample info with old and new names is available [here](http://149.156.177.112/projects/ifpan-chwastek-cells/analysis/full_sample_info.csv)

