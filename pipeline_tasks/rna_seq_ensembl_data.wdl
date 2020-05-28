workflow rna_seq_ensembl_data_workflow {

  meta {
    keywords: '{"keywords": ["ensembl dataset", "hisat2", "extract splice sites"]}'
    name: 'rna_seq_ensembl_data'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'Preparation of data from ensembl dataset'
    changes: '{"latest": "no changes"}'

    input_release_version: '{"name": "release_version", "type": "String", "description": "Ensembl release version"}'
    input_organism_name: '{"name": "organism_name", "type": "String", "description": "Name of the organism in Latin"}'

    output_ref_genome_index: '{"name": "ref_genome_index", "type": "Array[File]", "copy": "True", "description": "Reference genome indexed"}'
    output_gtf_file: '{"name": "gtf_file", "type": "File", "copy": "True", "description": "GTF file from ensembl dataset"}'
    output_splicesites_file: '{"name": "splicesites_file", "type": "File", "copy": "True", "description": "Splicesites file"}'
    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_log: '{"name": "Standard err", "tye": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'

  }

  call rna_seq_ensembl_data

}

task rna_seq_ensembl_data {

  String chromosome_name = "primary_assembly"
  String release_version
  String organism_name

  String genome_basename = "no_basename"
  String gtf_basename = sub(organism_name, " ", "_") + "_gtf"

  String task_name = "rna_seq_ensembl_data"
  String task_version = "latest"
  Int? index
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String docker_image = "intelliseqngs/hisat2:1.0.0"

  command <<<

  task_name="${task_name}"; task_name_with_index="${task_name_with_index}"; task_version="${task_version}"; task_docker="${docker_image}"
  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/after-start.sh)

    ORGANISM_NAME=$(echo ${organism_name} | sed -e 's/ /_/g' | awk '{print tolower($0)}')
    GENOME_FTP_PATH="ftp://ftp.ensembl.org/pub/release-${release_version}/fasta/$ORGANISM_NAME/dna/*.dna_sm.${chromosome_name}.fa.gz"
    GENOME_FTP_PATH_TOPLEVEL=$(echo $GENOME_FTP_PATH | sed -e 's/primary_assembly/toplevel/g')

    # check if ftp path exists (if exists exit status equal to 8, if does not exists exit status equal to 0)
    wget -q --spider $GENOME_FTP_PATH
    STATUS=$?
    if [ $STATUS -ne 0 ]
    then
        wget -O - $GENOME_FTP_PATH | gunzip > ${genome_basename}.fa
    else
        wget -O - $GENOME_FTP_PATH_TOPLEVEL | gunzip > ${genome_basename}.fa
    fi

    set -e pipefail

    hisat2-build ${genome_basename}.fa ${genome_basename}
    mkdir genome
    mv ${genome_basename}.* genome

    GTF_FTP_PATH="ftp://ftp.ensembl.org/pub/release-${release_version}/gtf/$ORGANISM_NAME/*.${release_version}.gtf.gz"

    wget -O - $GTF_FTP_PATH | gunzip > ${gtf_basename}.gtf
    extract_splice_sites.py ${gtf_basename}.gtf > ${gtf_basename}.splicesites

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/before-finish.sh)
  >>>

  runtime {

    docker: docker_image
    memory: "1G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    Array[File] ref_genome_index = glob("genome/*")
    File gtf_file = "${gtf_basename}.gtf"
    File splicesites_file = "${gtf_basename}.splicesites"
    String genome_index = "${genome_basename}"

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

  }

}
