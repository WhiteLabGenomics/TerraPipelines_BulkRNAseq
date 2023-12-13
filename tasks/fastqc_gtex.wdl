version 1.0
task fastqc {
    input {
        File fastq1
        File fastq2

        Float memory = 4
        Int disk_space = 5*ceil(size(fastq1, "GiB")) + 8
        Int num_threads = 1
        Int num_preempt = 0

        String fastq1_name = sub(sub(basename(fastq1), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )
        String fastq2_name = sub(sub(basename(fastq2), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )
    }

    command <<<
        set -euo pipefail
        fastqc ${fastq1} ${fastq2} \
            --threads ${num_threads} \
            --outdir .
        unzip -p ${fastq1_name}_fastqc.zip ${fastq1_name}_fastqc/fastqc_data.txt | gzip > ${fastq1_name}.fastqc_data.txt.gz
        unzip -p ${fastq2_name}_fastqc.zip ${fastq2_name}_fastqc/fastqc_data.txt | gzip > ${fastq2_name}.fastqc_data.txt.gz
    >>>

    output {
        File fastq1_fastqc_html = "${fastq1_name}_fastqc.html"
        File fastq1_fastqc_zip =  "${fastq1_name}_fastqc.zip"
        File fastq1_fastqc_data = "${fastq1_name}.fastqc_data.txt.gz"
        File fastq2_fastqc_html = "${fastq2_name}_fastqc.html"
        File fastq2_fastqc_zip =  "${fastq2_name}_fastqc.zip"
        File fastq2_fastqc_data = "${fastq2_name}.fastqc_data.txt.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow fastqc_workflow {
    input {
        File fastq1
        File fastq2
    }

    call fastqc {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2
    }

    output {
        File htmlReport_fastq1 = fastqc.fastq1_fastqc_html
        File reportZi_fastq1 = fastqc.fastq1_fastqc_zip
        File fastqc_data_fastq1 = fastqc.fastq1_fastqc_data
        File htmlReport_fastq2 = fastqc.fastq2_fastqc_html
        File reportZi_fastq2 = fastqc.fastq2_fastqc_zip
        File fastqc_data_fastq2 = fastqc.fastq2_fastqc_data
    }
}

