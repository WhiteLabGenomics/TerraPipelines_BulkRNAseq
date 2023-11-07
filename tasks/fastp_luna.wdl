version 1.0
# Adapter clipping
# https://github.com/OpenGene/fastp
workflow run_Fastp {
  call Fastp
}

task Fastp {
  input {
    File fastq1
    File fastq2
    String output_prefix
    File adapter_fasta = "gs://gcp-public-data--broad-references/RNA/resources/Illumina_adapters.fasta"

    String docker = "us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500"
    Int memory_mb =  ceil(1.5*size(fastq1, "MiB")) + 8192 # Experimentally determined formula for memory allocation
    Int disk_size_gb = 5*ceil(size(fastq1, "GiB")) + 128
    File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
    Int cpu=4
  }

  command {
    bash ${monitoring_script} > monitoring.log &

    fastp --in1 ${fastq1} --in2 ${fastq2} --out1 ${output_prefix}_read1.fastq.gz --out2 ${output_prefix}_read2.fastq.gz \
    --json ${output_prefix}_fastp.json \
    --html ${output_prefix}_fastp.html \
    --disable_quality_filtering \
    --disable_length_filtering \
    --adapter_fasta ${adapter_fasta} \
    --thread ${cpu}
  }
  

  runtime {
    docker: docker
    memory: "~{memory_mb} MiB"
    cpu: cpu
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: 0
    maxRetries: 2
  }

  output {
    File monitoring_log = "monitoring.log"
    File fastq1_clipped = output_prefix + "_read1.fastq.gz"
    File fastq2_clipped = output_prefix + "_read2.fastq.gz"
    File raport_json = output_prefix + "_fastp.json"
    File raport_html = output_prefix + "_fastp.html"
  }
}

workflow fastp_workflow {
    input {
        File fastq1
        File fastq2
        String output_prefix
    }

    call Fastp {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            output_prefix = output_prefix
    }

    output {
    File monitoring_log = Fastp.monitoring_log
    File fastq1_clipped = Fastp.fastq1_clipped
    File fastq2_clipped = Fastp.fastq2_clipped
    File raport_json = Fastp.raport_json
    File raport_html = Fastp.raport_html
    }
}

