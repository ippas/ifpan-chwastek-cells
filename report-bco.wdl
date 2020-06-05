workflow report_bco_workflow {


  meta {
    keywords: '{"keywords": []}'
    name: 'report_bco'
    author: 'https://gitlab.com/moni.krzyz'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'report_bco \n Generic text for task'
    changes: '{"latest": "no changes"}'
    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_err: '{"name": "Standard err", "type": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'
  }

  call report_bco

}

task report_bco {

  String task_name = "report_bco"
  String task_version = "latest"
  String docker_image = "intelliseqngs/reports:2.3.0"
  File bco_json
  String sample_id = "no_id_provided"

  command <<<
  task_name="${task_name}"; task_version="${task_version}"; task_docker="${docker_image}"
  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.0/after-start.sh)

  /intelliseqtools/generate-report.sh --json bco=${bco_json} --template /intelliseqtools/templates/bco-v1/content.jinja --name ${sample_id}.bco_report

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.0/before-finish.sh)
  >>>

  runtime {

    docker: docker_image
    memory: "1G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

    File bco_report_pdf = "${sample_id}.bco_report.pdf"
    File bco_report_odt = "${sample_id}.bco_report.odt"
    File bco_report_docx = "${sample_id}.bco_report.docx"
    File bco_report_html = "${sample_id}.bco_report.html"

  }
}
