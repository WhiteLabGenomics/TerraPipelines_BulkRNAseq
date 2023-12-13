version 1.0

import "../tasks/fastqc.wdl" as fastqc_task
import "../tasks/fastp.wdl" as fastp_task
import "../tasks/star.wdl" as star_task
import "../tasks/rnaseqc2.wdl" as rnaseqc2_task
import "../tasks/rsem.wdl" as rsem_task


workflow RNA_preprocessing_pipeline {

  input {

    File fastq1
    File fastq2
    String sample_id

    #annotation GTF
    File genes_gtf="gs://ccle_default_params/references_gtex_gencode.v29.GRCh38.ERCC.genes.collapsed_only.gtf"

    #star_task index
    File star_index="gs://ccle_default_params/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh100.tar.gz"

    #rsem index
    File rsem_reference="gs://ccle_default_params/rsem_reference_GRCh38_gencode29_ercc.tar.gz"

  }

  call fastp_task.Fastp as fastp {
      input:
        fastq1 = fastq1,
        fastq2 = fastq2,
        output_prefix = sample_id
    }

  call fastqc_task.fastqc as cleaned_fastqc1 {
    input:
      seqFile=fastp.fastq1_clipped,
  }

  call fastqc_task.fastqc as cleaned_fastqc2 {
    input:
      seqFile=fastp.fastq2_clipped,
  }

  call star_task.star as star {
    input:
      prefix=sample_id,
      fastq1=fastp.fastq1_clipped,
      fastq2=fastp.fastq2_clipped,
      star_index=star_index
  }

  call rnaseqc2_task.rnaseqc2 as rnaseqc2 {
    input:
      bam_file=star.bam_file,
      genes_gtf=genes_gtf,
      sample_id=sample_id
  }

  call rsem_task.rsem as rsem {
    input:
      transcriptome_bam=star.transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference,
  }


  output {
    #fastp
    File monitoring_log = fastp.monitoring_log
    File fastq1_clipped = fastp.fastq1_clipped
    File fastq2_clipped = fastp.fastq2_clipped
    File raport_json = fastp.raport_json
    File raport_html = fastp.raport_html

    #fastqc cleaned data
    File cleaned_htmlReport1 = cleaned_fastqc1.htmlReport
    File cleaned_reportZip1 = cleaned_fastqc1.reportZip
    File cleaned_htmlReport2 = cleaned_fastqc2.htmlReport
    File cleaned_reportZip2 = cleaned_fastqc2.reportZip

    #star
    File bam_file=star.bam_file
    File bam_index=star.bam_index
    File transcriptome_bam=star.transcriptome_bam
    File chimeric_junctions=star.chimeric_junctions
    File chimeric_bam_file=star.chimeric_bam_file
    File read_counts=star.read_counts
    File junctions=star.junctions
    File junctions_pass1=star.junctions_pass1
    Array[File] logs=star.logs

    #rnaseqc
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr

    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms
  }

}

