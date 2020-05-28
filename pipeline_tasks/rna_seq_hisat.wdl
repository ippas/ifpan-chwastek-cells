workflow rna_seq_hisat_workflow {

  meta {
    keywords: '{"keywords": ["alignment"]}'
    name: 'rna_seq_hisat'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'Alignment with hisat2'
    changes: '{"latest": "no changes"}'

    input_fastq_1: '{"name": "fastq_1", "type": "File", "description": "Fastq file 1 (paired-end)"}'
    input_fastq_2: '{"name": "fastq_2", "type": "File", "description": "Fastq file 1 (paired-end)"}'
    input_ref_genome_index: '{"name": "ref_genome_index", "type": "Array[File]", "description": "reference genome indexed"}'
    input_splicesites_file: '{"name": "splicesites_file", "type": "File", "description": "splicesites file"}'
    input_sample_id: '{"name": "sample id", "type": "String", "description": "identifier of sample"}'

    output_bam_file: '{"name": "bam_file", "type": "File", "copy": "True", "description": "bam file"}'
    output_bam_bai_file: '{"name": "bam_bai_file", "type": "File", "copy": "True", "description": "bam bai file"}'
    output_summary: '{"name": "summary", "type": "File", "copy": "True", "description": "summary file"}'
    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_log: '{"name": "Standard err", "type": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'

  }

  call rna_seq_hisat

}

task rna_seq_hisat {

  File fastq_1
  File fastq_2
  Array[File] ref_genome_index
  File splicesites_file
  String genome_basename = "no_basename"

  String num_cpu = "4"
  String sample_id = "no_id_provided"

  String task_name = "rna_seq_hisat"
  String task_version = "latest"
  Int? index
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String docker_image = "intelliseqngs/hisat2:1.0.0"

  command <<<
  task_name="${task_name}"; task_name_with_index="${task_name_with_index}"; task_version="${task_version}"; task_docker="${docker_image}"
  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/after-start.sh)

    set -e pipefail

    mkdir ${genome_basename}
    ln -s ${sep=" " ref_genome_index} ${genome_basename}

    hisat2 -t --known-splicesite-infile ${splicesites_file} --dta-cufflinks -x ${genome_basename}/${genome_basename} \
    --rna-strandness FR -1 ${fastq_1} -2 ${fastq_2} -p ${num_cpu} --summary-file ${sample_id}.txt \
    | samtools sort -@ ${num_cpu} -O BAM -o ${sample_id}.bam
    samtools index ${sample_id}.bam ${sample_id}.bam.bai

    rm ${fastq_1} ${fastq_2} ${splicesites_file} ${sep=" " ref_genome_index}

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/before-finish.sh)
  >>>

  runtime {

    docker: docker_image
    memory: "1G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    File bam_file = "${sample_id}.bam"
    File bam_bai_file = "${sample_id}.bam.bai"
    File summary = "${sample_id}.txt"

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

  }

}
