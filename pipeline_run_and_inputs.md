##### Analysis steps:

1. Preprocessing was performed with [intelliseq workflow](rna-seq-paired-end.wdl). Example bco object produced by this pipeline is available [here](bco.json)

2. Inputs were prepared with the following bash script and processed in batches of 10-14 (due to RAM restrictions):

```
#! /bin/bash
echo {\"rna_seq_paired_end.fastqs\": [
ls */*/*1.fq.gz | xargs -i bash -c 'echo \" {} \"',
ls */*/*2.fq.gz | xargs -i bash -c 'echo \" {} \"',
echo ], \"rna_seq_paired_end.organism_name\": \"Homo Sapiens\", \"rna_seq_paired_end.release_version\": \"100\" }
```
*note, this requires some manual postprocessing: delete excess commas and all spaces deleted: 
`sed -i 's/ //g' input.json`, but be carefull to keep the space in Homo sapiens*

3. Command to run the whole pipeline (in a screen session):
`java -jar /opt/tools/cromwell/cromwell-44.jar run rna_seq_paired.wdl -i input_x.json > rna_paired_log_x`

4. Cufflinks abundancies files from all batches were collected and normalised with cuffnorm together in [this container](link to container)
