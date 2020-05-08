## Osteoarthritis model RNA-Seq 

Patients cells treated with 3 compounds (4 h, 10 ng/ml), 2-TNFa, 3-LPS, 4-IFg.

* [list of samples](samples.csv)

### Analysis steps:

1. Preprocessing was performed with [intelliseq workflow](), whole pipeline with all of the tasks is avaiable in this repo.

2. Inputs were prepared with the following bash script:

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
which was manually post-processed (3 excess commas deleted - these at the end of each file list)


3. Command to run the whole pipeline:

############`less samples.txt | xargs -i bash -c 'java -jar /opt/tools/cromwell/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/qc-fq-fastqc/latest/qc-fq-fastqc.wdl -i {}-input.json > {}-log.txt' &` *the '&' sign at the end of line tells bash to run whatever command in the background

