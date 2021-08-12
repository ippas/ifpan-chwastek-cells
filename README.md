## RNA-Seq in Ostheoarthritis model

Patients cells treated with 3 compounds (4 h, 10 ng/ml), 1-ctrl 2-IFNg, 3-TNFa, 4-LPS

### Methods:

All samples were checked for quality with fastQC v0.11.9 and aligned to a GRCh38 *Homo sapiens* reference genome (release 100 from Ensembl database) with hisat2 2.2.0. Cufflinks v 2.2.1 package and GTF from the Ensembl gene database were used to quantify (cuffquant) and normalize (cuffnorm) transcripts to fpkms (Fragments Per Kilobase of transcript per Million fragments mapped). All statisical analyses were performed with R software v3.6.


*note: a few samples were mislabeled in this experiment. The labels were corrected by switching sample names in sample info and are available in results folder
