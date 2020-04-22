## Osteoarthritis model RNA-Seq 

Patients cells treated with 3 compounds (4 h, 10 ng/ml), 2-TNFa, 3-LPS, 4-IFg.

* [list of samples](samples.csv)

### Analysis steps:

1. Sample quality control:

Fastqc was run with the [Intelliseq workflow](https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl) in cromwell/wdl in batches of 10 files. With the following code:

This code is to run intelliseq's task of fastqc:

This line lists all sample names based on the file_1.fq.gz file_2.fq.gz naming convention:

`ls */*/*1.fq.gz | xargs -i bash -c 'BASENAME2=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME2' > samples.txt`


To generate json input files (generated earlier with WOMtool) from names:
`less samples.txt | xargs -i bash -c 'echo "{\"qc_fq_fastqc_workflow.qc_fq_fastqc.fastq_1\":\"{}_1.fq.gz\",\"qc_fq_fastqc_workflow.qc_fq_fastqc.fastq_2\":\"{}_2.fq.gz\"}">{}-input.json'`


To run the workflow on an example file: 
`java -jar /opt/tools/cromwell/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/qc-fq-fastqc/latest/qc-fq-fastqc.wdl -i ./fq/H011/H011-input.json > H011-log.json`


To run the workflow on all files: `less samples.txt | xargs -i bash -c 'java -jar /opt/tools/cromwell/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/qc-fq-fastqc/latest/qc-fq-fastqc.wdl -i {}-input.json > {}-log.txt' &` *the '&' sign at the end of line tells bash to run whatever command in the background


* alternatively I am trying to use the whole pipeline

Multiqc was used to generate the final report, in this [docker container](https://hub.docker.com/r/ewels/multiqc). Command: `docker run --rm -v $PWD:/data ewels/multiqc:latest multiqc /data -o /data`.

+ add link to fastqc report
