version 1.0

import "star.wdl" as star_v1
import "rnaseqc2.wdl" as rnaseqc2_v1
import "rsem.wdl" as rsem_v1
import "star_fusion.wdl" as star_fusion
import "rnaseq-germline-snps-indels.wdl" as rnaseq_mutations


workflow RNA_pipeline {

  input {

    File fastq1
    File fastq2
    String sample_id

    # mutations
    File refFasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File refFastaIndex="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File refDict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

    Array[File] knownVcfs=["gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
                            "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"]
    Array[File] knownVcfsIndices=['gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi',
                                  "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"]
    File dbSnpVcf="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    File dbSnpVcfIndex="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
    File annotationsGTF="gs://gcp-public-data--broad-references/hg38/v0/gencode.v27.primary_assembly.annotation.gtf"#gs://gatk-test-data/intervals/star.gencode.v19.transcripts.patched_contigs.gtf
    # TODO: create a test.wdl.json
    # TODO: create a local_test.wdl.json
    # TODO: add in other side experimental genomes that could be in our files (like ERCC)
    # TODO: put all at the same gencode version (latest)
    #star_v1
    File star_index= "gs://ccle_default_params/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v29_oh100.tar.gz"

    #star fusion
    Array[File] ctat_genome_lib_build_dir_files=[
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/AnnotFilterRule.pm",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/blast_pairs.idx",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/blast_pairs.idx.prev.1553723931",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/fusion_annot_lib.idx",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/pfam_domains.dbm",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.cds",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.mini.sortu",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.pep",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.prot_info.dbm",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.fai",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm"
    ]
    Array[File] ref_genome_fa_star_idx_files=[
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/Genome",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SA",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SAindex",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/build.ok",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrLength.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrName.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrStart.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/exonGeTrInfo.tab",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/exonInfo.tab",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/geneInfo.tab",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbInfo.txt",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab",
      "gs://ccle_default_params/references/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/transcriptInfo.tab"
    ]

    #rsem
    File rsem_reference="gs://ccle_default_params/rsem_reference_GRCh38_gencode29_ercc.tar.gz"
    
    #rnaseqc2_v1
    File genes_gtf="gs://ccle_default_params/references_gtex_gencode.v29.GRCh38.ERCC.genes.collapsed_only.gtf"

    #rna_mutect2
    #Boolean run_funcotator=false
    #String gatk_docker="broadinstitute/gatk:4.2.2.0"
    #Array[File] knownVcfs
    #Array[File] knownVcfsIndices
    #File dbSnpVcf
    #File dbSnpVcfIndex
    #File ref_fasta
    #File ref_fai
    #File ref_dict
    #String source
    #Boolean compress_vcfs=true
    #Boolean filter_funcotations=false
    #Boolean funco_compress=true
    #String funco_output_format="VCF"
    #String funco_reference_version="hg38"
    #Boolean funco_use_gnomad_AF=true
    #Int scatter_count=20
    #String gcs_project_for_requester_pays="broad-firecloud-ccle"
    #File intervals="gs://gcp-public-data--broad-references/hg38/v0/wgs_coverage_regions.hg38.interval_list"
    #File gnomad = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
    #File gnomad_idx = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
    #String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
    #Boolean make_bamout=false
    #File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
    #File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
  }


  call star_v1.star as star {
    input:
      prefix=sample_id,
      fastq1=fastq1,
      fastq2=fastq2,
      star_index=star_index
  }

  call rnaseqc2_v1.rnaseqc2 as rnaseqc2 {
    input:
      bam_file=star.bam_file,
      genes_gtf=genes_gtf,
      sample_id=sample_id
  }

  call rnaseq_mutations.RNAseq as rnaseq_mutations {
    input:
      inputBam=star.bam_file,
      sampleName=sample_id,
      refFasta=refFasta,
      refFastaIndex=refFastaIndex,
      refDict=refDict,
      knownVcfs=knownVcfs,
      knownVcfsIndices=knownVcfsIndices,
      dbSnpVcf=dbSnpVcf,
      dbSnpVcfIndex=dbSnpVcfIndex,
      annotationsGTF=annotationsGTF
  }
  # TODO: annotate the mutations
  # bcftools? (need to see some vcfs)
  # cnn_filter (once trainned on )
  # TODO: train CNN filter from GATK on RNAseq data
  # TODO: train a NN (enformer-like) on predicting batch informations (sequencer, reagents, RNAtype, species, iscancer, ERCC spike in) directly from fastq data
  # remove filtered
  # open cravat (need to update the annotators)
  # 
  # build a tool to extract metadata from the QC files and bams
  # build a tool to merge over multiple samples given a list of sample names
  # get the tool to move data in Terra from the terra bucket to our bucket
  # add a cleanup step for the terra jobs. 

  call rsem_v1.rsem as rsem {
    input:
      transcriptome_bam=star.transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference,
      is_stranded="false",
      paired_end="true"
  }

  # TODO: give junctions directly to star fusion
  call star_fusion.StarFusion as StarFusion {
    input:
      left_fastq=fastq1,
      right_fastq=fastq2,
      prefix=sample_id,
      ctat_genome_lib_build_dir_files=ctat_genome_lib_build_dir_files,
      ref_genome_fa_star_idx_files=ref_genome_fa_star_idx_files
  }

  output {
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
    #StarFusion
    File fusion_predictions=StarFusion.fusion_predictions
    File fusion_predictions_abridged=StarFusion.fusion_predictions_abridged
    # mutations
    File variant_filtered_vcf=rnaseq_mutations.variant_filtered_vcf
    File variant_filtered_vcf_index=rnaseq_mutations.variant_filtered_vcf_index

  }
}

