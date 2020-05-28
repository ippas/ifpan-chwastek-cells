import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-ensembl-data@1.0.3/src/main/wdl/tasks/rna-seq-ensembl-data/latest/rna-seq-ensembl-data.wdl" as rna_seq_ensembl_data_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-hisat@1.0.1/src/main/wdl/tasks/rna-seq-hisat/latest/rna-seq-hisat.wdl" as rna_seq_hisat_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-cuffquant@1.0.1/src/main/wdl/tasks/rna-seq-cuffquant/latest/rna-seq-cuffquant.wdl" as rna_seq_cuffquant_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-cuffnorm@1.0.1/src/main/wdl/tasks/rna-seq-cuffnorm/latest/rna-seq-cuffnorm.wdl" as rna_seq_cuffnorm_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-concat-summary@1.0.1/src/main/wdl/tasks/rna-seq-concat-summary/latest/rna-seq-concat-summary.wdl" as rna_seq_concat_summary_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-fastqc@1.0.1/src/main/wdl/tasks/rna-seq-fastqc/latest/rna-seq-fastqc.wdl" as rna_seq_fastqc_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-qc-stats@1.0.1/src/main/wdl/tasks/rna-seq-qc-stats/latest/rna-seq-qc-stats.wdl" as rna_seq_qc_stats_task
#import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-fastqc-overrep@1.0.1/src/main/wdl/tasks/rna-seq-fastqc-overrep/latest/rna-seq-fastqc-overrep.wdl" as rna_seq_fastqc_overrep_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-merge-bco@1.0.0/src/main/wdl/tasks/rna-seq-merge-bco/latest/rna-seq-merge-bco.wdl" as rna_seq_merge_bco_task
#import "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/report-bco/latest/report-bco.wdl" as report_bco_task

workflow rna_seq_paired_end {

  meta {
    name: 'rna_seq_paired_end'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: '## RNA-Seq (paired-end)'
    changes: '{"latest": "no changes"}'

    input_fastqs_1: '{"name": "fastq 1", "type": "Array[File]", "description": "first fastq file"}'
    input_fastqs_2: '{"name": "fastq 2", "type": "Array[File]", "description": "second fastq file"}'
    input_samples_ids: '{"name": "samples_ids", "type": "Array[String]", "description": "identifiers of samples"}'
    input_organism_name: '{"name": "organism_name",  "type": "String", "description": "name of the organism in Latin"}'
    input_release_version: '{"name": "release_version",  "type": "String", "description": "ensembl release version"}'

    output_fastqc_zip_files: '{"name": "fastqc_zip_files", "type": "Array[Array[File]]", "copy": "True", "description": "zip files from FastQC"}'
    output_basic_statistics_excel: '{"name": "basic_statistics_excel", "type": "File", "copy": "true", "description": "basic statistics in Excel", "constraints": {"extension": ["xlsx"]}}'
    output_basic_statistics_csv: '{"name": "basic_statistics_csv", "type": "File", "copy": "true", "description": "basic statistics in CSV", "constraints": {"extension": ["csv"]}}'
    output_overrepresented: '{"name": "overrepresented", "type": "File", "copy": "true", "description": "overrepresented in Excel", "constraints": {"extension": ["xlsx"]}}'
    output_overrepresented_csv: '{"name": "overrepresented_csv", "type": "File", "copy": "true", "description": "overrepresented in CSV", "constraints": {"extension": ["csv"]}}'
    output_bam_file: '{"name": "bam_file", "type": "Array[File]", "copy": "True", "description": "bam files"}'
    output_bam_bai_file: '{"name": "bam_bai_file", "type": "Array[File]", "copy": "True", "description": "bam bai files"}'
    output_summary_excel_file: '{"name": "summary_excel_file", "type": "File", "copy": "true", "description": "summary file in Excel", "constraints": {"extension": ["xlsx"]}}'
    output_summary_csv_file: '{"name": "summary_csv_file", "type": "File", "copy": "true", "description": "summary file in CSV", "constraints": {"extension": ["csv"]}}'
    output_abundances_file: '{"name": "abundances_file", "type": "Array[File]", "copy": "True", "description": "abundances files"}'
    output_cuffnorm_output: '{"name": "cuffnorm_output", "type": "Array[File]", "copy": "True", "description": "cuffnorm outputs"}'
  }

  Array[File] fastqs_1
  Array[File] fastqs_2
  Array[String] samples_ids
  String organism_name
  String release_version
  String genome_basename = sub(organism_name, " ", "_") + "_genome"
  String chromosome_name = "primary_assembly"
  String summary_file_name = "summary"
  String pipeline_name = "rna_seq_paired_end"
  String pipeline_version = "latest"

  scatter (index in range(length(fastqs_1))) {
    call rna_seq_fastqc_task.rna_seq_fastqc {
        input:
            fastq_1 = fastqs_1[index],
            fastq_2 = fastqs_2[index],
            sample_id = samples_ids[index],
            index = index
    }
  }

  call rna_seq_qc_stats_task.rna_seq_qc_stats {
    input:
        zip_files = rna_seq_fastqc.zip_files
  }

#  call rna_seq_fastqc_overrep_task.rna_seq_fastqc_overrep {
#    input:
#        zip_files = rna_seq_fastqc.zip_files
#  }

  call rna_seq_ensembl_data_task.rna_seq_ensembl_data {
    input:
        release_version = release_version,
        chromosome_name = chromosome_name,
        organism_name = organism_name,
        genome_basename = genome_basename
  }

  scatter (index in range(length(fastqs_1))) {
    call rna_seq_hisat_task.rna_seq_hisat {
        input:
            fastq_1 = fastqs_1[index],
            fastq_2 = fastqs_2[index],
            sample_id = samples_ids[index],
            splicesites_file = rna_seq_ensembl_data.splicesites_file,
            ref_genome_index = rna_seq_ensembl_data.ref_genome_index,
            genome_basename = genome_basename,
            index = index
    }
  }

  call rna_seq_concat_summary_task.rna_seq_concat_summary {
    input:
        summary_files = rna_seq_hisat.summary,
        summary_file_name = summary_file_name
  }

  scatter (index in range(length(fastqs_1))) {
    call rna_seq_cuffquant_task.rna_seq_cuffquant {
        input:
            gtf_file = rna_seq_ensembl_data.gtf_file,
            bam_file = rna_seq_hisat.bam_file[index],
            sample_id = samples_ids[index],
            index = index
    }
  }

  call rna_seq_cuffnorm_task.rna_seq_cuffnorm {
    input:
        gtf_file = rna_seq_ensembl_data.gtf_file,
        abundances_files = rna_seq_cuffquant.abundances_file
  }

# Merge bco, stdout, stderr files
  Array[File] bco_tasks = [rna_seq_qc_stats.bco, rna_seq_ensembl_data.bco, rna_seq_concat_summary.bco, rna_seq_cuffnorm.bco]
  Array[File] stdout_tasks = [rna_seq_qc_stats.stdout_log, rna_seq_ensembl_data.stdout_log, rna_seq_concat_summary.stdout_log, rna_seq_cuffnorm.stdout_log]
  Array[File] stderr_tasks = [rna_seq_qc_stats.stderr_log, rna_seq_ensembl_data.stderr_log, rna_seq_concat_summary.stderr_log, rna_seq_cuffnorm.stderr_log]

  Array[Array[File]] bco_scatters = [bco_tasks, rna_seq_fastqc.bco, rna_seq_hisat.bco, rna_seq_cuffquant.bco]
  Array[Array[File]] stdout_scatters = [stdout_tasks, rna_seq_fastqc.stdout_log, rna_seq_hisat.stdout_log, rna_seq_cuffquant.stdout_log]
  Array[Array[File]] stderr_scatters = [stderr_tasks, rna_seq_fastqc.stderr_log, rna_seq_hisat.stderr_log, rna_seq_cuffquant.stderr_log]

  Array[File] bco_array = flatten(bco_scatters)
  Array[File] stdout_array = flatten(stdout_scatters)
  Array[File] stderr_array = flatten(stderr_scatters)

  call rna_seq_merge_bco_task.rna_seq_merge_bco {
    input:
        bco_array = bco_array,
        stdout_array = stdout_array,
        stderr_array = stderr_array,
        pipeline_name = pipeline_name,
        pipeline_version = pipeline_version
  }

  output {
    Array[Array[File]] fastqc_zip_files = rna_seq_fastqc.zip_files
    File basic_statistics_excel = rna_seq_qc_stats.basic_statistics_excel
    File basic_statistics_csv = rna_seq_qc_stats.basic_statistics_csv
#    File overrepresetned = rna_seq_fastqc_overrep.overrepresented
#    File overrepresetned_csv = rna_seq_fastqc_overrep.overrepresented_csv
    Array[File] bam_file = rna_seq_hisat.bam_file
    Array[File] bam_bai_file = rna_seq_hisat.bam_bai_file
    File summary_excel_file = rna_seq_concat_summary.summary_excel_file
    File summary_csv_file = rna_seq_concat_summary.summary_csv_file
    Array[File] abundances_file = rna_seq_cuffquant.abundances_file
    Array[File] cuffnorm_output = rna_seq_cuffnorm.cuffnorm_output
  }
}
