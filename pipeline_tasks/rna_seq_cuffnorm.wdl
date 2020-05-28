workflow rna_seq_cuffnorm_workflow {

  meta {
    keywords: '{"keywords": ["cuffnorm"]}'
    name: 'rna_seq_cuffnorm'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'Generic text for task'
    changes: '{"latest": "no changes"}'

    input_sample_id: '{"name": "sample id", "type": "String", "description": "identifier of sample"}'
    input_gtf_file: '{"name": "gtf_file", "type": "File", "description": "gtf file"}'
    input_abundances_files: '{"name": "bam_file", "type": "Array[File]", "description": "abundances file"}'

    output_cuffnorm_output: '{"name": "cuffnorm_output", "type": "File", "copy": "True", "description": "cuffnorm output"}'
    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_log: '{"name": "Standard err", "type": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'

  }

  call rna_seq_cuffnorm

}

task rna_seq_cuffnorm {

  File gtf_file
  Array[File] abundances_files

  String num_cpu = "4"
  String sample_id = "no_id_provided"

  String task_name = "rna_seq_cuffnorm"
  String task_version = "latest"
  Int? index
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String docker_image = "intelliseqngs/cufflinks:1.0.0"

  command <<<
  task_name="${task_name}"; task_name_with_index="${task_name_with_index}"; task_version="${task_version}"; task_docker="${docker_image}"
  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/after-start.sh)

  set -e pipefail

  ### Cuffnorm
  cuffnorm -p ${num_cpu} -v -library-norm-method classic-fpkm -o cuffnorm \
        ${gtf_file} ${sep=" " abundances_files}

  python <<CODE

  import pandas as pd
  import os

  sample_names = pd.read_csv("cuffnorm/samples.table", sep="\t")
  count_and_fpkm_files = [f for f in os.listdir("cuffnorm/") if f.endswith(('.count_table', '.fpkm_table'))]
  os.chdir("cuffnorm/")

  for i in range(len(count_and_fpkm_files)):
      data = pd.read_csv(count_and_fpkm_files[i], sep="\t", index_col=0)
      if data.empty:
          continue
      columns_name = []
      for j in range(data.shape[1]):
          columns_name.append(sample_names.iloc[j,1].split("/")[-1].split(".cxb")[0])
      data.columns = columns_name
      data.to_csv(count_and_fpkm_files[i], sep="\t")

  CODE

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/before-finish.sh)

  >>>

  runtime {

    docker: docker_image
    memory: "1G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    Array[File] cuffnorm_output = glob("cuffnorm/*")

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

  }

}
