workflow rna_seq_cuffquant_workflow {

  meta {
    keywords: '{"keywords": ["cuffquant"]}'
    name: 'rna_seq_cuffquant'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'Generic text for task'
    changes: '{"latest": "no changes"}'

    input_sample_id: '{"name": "sample id", "type": "String", "description": "identifier of sample"}'
    input_gtf_file: '{"name": "gtf_file", "type": "File", "description": "gtf file"}'
    input_bam_file: '{"name": "bam_file", "type": "File", "description": "bam_file"}'

    output_abundances_file: '{"name": "abundances_file", "type": "File", "copy": "True", "description": "abundances file"}'
    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_log: '{"name": "Standard err", "type": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'

  }

  call rna_seq_cuffquant

}

task rna_seq_cuffquant {

  File gtf_file
  File bam_file

  String sample_id = "no_id_provided"
  String num_cpu = "4"

  String task_name = "rna_seq_cuffquant"
  String task_version = "latest"
  Int? index
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String docker_image = "intelliseqngs/cufflinks:1.0.0"

  command <<<
  task_name="${task_name}"; task_name_with_index="${task_name_with_index}"; task_version="${task_version}"; task_docker="${docker_image}"
  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/after-start.sh)

    set -e pipefail

    cuffquant -v -p ${num_cpu} --no-effective-length-correction -o cuffquant \
     ${gtf_file} ${bam_file}

    mv cuffquant/abundances.cxb cuffquant/${sample_id}.cxb

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/before-finish.sh)
  >>>

  runtime {

    docker: docker_image
    memory: "1G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    File abundances_file = "cuffquant/${sample_id}.cxb"

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

  }

}
