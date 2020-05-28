## Osteoarthritis model RNA-Seq 

Patients cells treated with 3 compounds (4 h, 10 ng/ml), 2-TNFa, 3-LPS, 4-IFg.

* [list of samples](samples.csv)

### Methods:

************** DO UZUPEŁNIENIA ********************


##### Analysis steps:

1. Preprocessing was performed with [intelliseq workflow](rna-seq-paired-end.wdl). Example bco object produced by this pipeline is available [here](bco.json)

2. Inputs were prepared with the following bash script and processed in batches of 10:

```
#! /bin/bash
echo {\"rna_seq_paired_end.samples_ids\":[
ls */*/*1.fq.gz | xargs -i bash -c 'BASENAME2=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1 | cut -d "/" -f 2); echo \"$BASENAME2\",'
echo ], \"rna_seq_paired_end.fastqs_2\": [
ls */*/*1.fq.gz | xargs -i bash -c 'echo \" {} \"',
echo ], \"rna_seq_paired_end.organism_name\": \"Homo Sapiens\", \"rna_seq_paired_end.fastqs_1\":[
ls */*/*2.fq.gz | xargs -i bash -c 'echo \" {} \"',
echo ], \"rna_seq_paired_end.release_version\": \"100\" }
```
which was manually post-processed (3 excess commas deleted - these at the end of each file list and all spaces deleted: 
`sed -i 's/ //g' input.json`, but be carefull to keep the space in Homo sapiens)


3. Command to run the whole pipeline (in a screen session):
`java -jar /opt/tools/cromwell/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/-/raw/dev/src/main/wdl/pipelines/rna-seq-paired-end/latest/rna-seq-paired-end.wdl -i input_x.json > rna_paired_log_x`
