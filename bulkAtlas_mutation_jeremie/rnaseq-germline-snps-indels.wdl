version 1.0

## Copyright Broad Institute, 2019
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v4) and related tools
##
## Requirements/expectations :
## - BAM 
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs. 

workflow RNAseq {

    input {
        File inputBam
        String sampleName

        File refFasta
        File refFastaIndex
        File refDict

        String gatk4_docker = "broadinstitute/gatk:latest"
        String gatk_path = "/gatk/gatk"
        String star_docker = "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"

        Array[File] knownVcfs
        Array[File] knownVcfsIndices

        File dbSnpVcf
        File dbSnpVcfIndex

        Int? minConfidenceForVariantCalling

        ## Inputs for STAR
        Int? readLength
        File? zippedStarReferences
        File annotationsGTF

        ## Optional user optimizations
        Int haplotypeScatterCount=6

        Int preemptible_tries=3
        # funcotator
        String? sequencing_center
        String? sequence_source
        String? funco_reference_version
        String funco_output_format="VCF"
        Boolean funco_compress=true
        Boolean funco_use_gnomad_AF=false
        File? funco_data_sources_tar_gz
        String? funco_transcript_selection_mode
        File? funco_transcript_selection_list
        Array[String]? funco_annotation_defaults
        Array[String]? funco_annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean funco_filter_funcotations=true
        String? funcotator_extra_args
    }

    Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
    
    call gtfToCallingIntervals {
        input:
            gtf = annotationsGTF,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            gatk_path = gatk_path,
            docker = gatk4_docker
    }

    call MarkDuplicates {
        input:
            input_bam = inputBam,
            base_name = sampleName + ".dedupped",
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    call SplitNCigarReads {
        input:
            input_bam = MarkDuplicates.output_bam,
            input_bam_index = MarkDuplicates.output_bam_index,
            base_name = sampleName + ".split",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    call BaseRecalibrator {
        input:
            input_bam = SplitNCigarReads.output_bam,
            input_bam_index = SplitNCigarReads.output_bam_index,
            recal_output_file = sampleName + ".recal_data.csv",
            dbSNP_vcf = dbSnpVcf,
            dbSNP_vcf_index = dbSnpVcfIndex,
            known_indels_sites_VCFs = knownVcfs,
            known_indels_sites_indices = knownVcfsIndices,
            ref_dict = refDict,
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }

    call ApplyBQSR {
        input:
            input_bam =  SplitNCigarReads.output_bam,
            input_bam_index = SplitNCigarReads.output_bam_index,
            base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            recalibration_report = BaseRecalibrator.recalibration_report,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    call ScatterIntervalList {
        input:
            interval_list = gtfToCallingIntervals.interval_list,
            scatter_count = haplotypeScatterCount,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }


    scatter (interval in ScatterIntervalList.out) {
        call HaplotypeCaller {
            input:
                input_bam = ApplyBQSR.output_bam,
                input_bam_index = ApplyBQSR.output_bam_index,
                base_name = sampleName + ".hc",
                interval_list = interval,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                ref_dict = refDict,
                dbSNP_vcf = dbSnpVcf,
                dbSNP_vcf_index = dbSnpVcfIndex,
                stand_call_conf = minConfidenceForVariantCalling,
                preemptible_count = preemptible_tries,
                docker = gatk4_docker,
                gatk_path = gatk_path
        }

    }

    call MergeVCFs {
        input:
            input_vcfs = HaplotypeCaller.output_vcf,
            input_vcfs_indexes =  HaplotypeCaller.output_vcf_index,
            output_vcf_name = sampleName + ".g.vcf.gz",
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }
    
    call VariantFiltration {
        input:
            input_vcf = MergeVCFs.output_vcf,
            input_vcf_index = MergeVCFs.output_vcf_index,
            base_name = sampleName + ".variant_filtered.vcf.gz",
            ref_fasta = refFasta,
            ref_fasta_index = refFastaIndex,
            ref_dict = refDict,
            preemptible_count = preemptible_tries,
            docker = gatk4_docker,
            gatk_path = gatk_path
    }

    Int disk_pad = 15
    Float file_size_multiplier=2.25

    call Funcotate {
        input:
            ref_fasta = refFasta,
            ref_fai = refFastaIndex,
            ref_dict = refDict,
            input_vcf = VariantFiltration.output_vcf,
            input_vcf_idx = VariantFiltration.output_vcf_index,
            # could be none
            reference_version = select_first([funco_reference_version, "hg38"]),
            output_file_base_name = basename(VariantFiltration.output_vcf, ".vcf") + ".annotated",
            output_format = funco_output_format,
            compress = funco_compress,
            use_gnomad = funco_use_gnomad_AF,
            data_sources_tar_gz = funco_data_sources_tar_gz,

            sequencing_center = sequencing_center,
            sequence_source = sequence_source,
            transcript_selection_mode = funco_transcript_selection_mode,
            transcript_selection_list = funco_transcript_selection_list,
            annotation_defaults = funco_annotation_defaults,
            annotation_overrides = funco_annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            interval_list=gtfToCallingIntervals.interval_list, # maybe not the right one
            filter_funcotations = funco_filter_funcotations,
            extra_args = funcotator_extra_args,
            disk_space = ceil(size(VariantFiltration.output_vcf, "GB") * file_size_multiplier)  + funco_tar_size + disk_pad,
            gatk_docker=gatk4_docker,
            machine_mem=4000,
            preemptible=2,
            max_retries=2,
            cpu=2,

    }

    output {
        File recalibrated_bam = ApplyBQSR.output_bam
        File recalibrated_bam_index = ApplyBQSR.output_bam_index
        File merged_vcf = MergeVCFs.output_vcf
        File merged_vcf_index = MergeVCFs.output_vcf_index
        File variant_filtered_vcf = VariantFiltration.output_vcf
        File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
        File funcotated_output_file = Funcotate.funcotated_output_file
        File funcotated_output_file_index = Funcotate.funcotated_output_file_index
    }
}

######### TASKS #########
######### TASKS #########
######### TASKS #########
######### TASKS #########

task gtfToCallingIntervals {
    input {
        File gtf
        File ref_dict

        String output_name = basename(gtf, ".gtf") + ".exons.interval_list"

        String docker
        String gatk_path
        Int preemptible_count
    }



    command <<<
        set -e

        Rscript --no-save -<<'RCODE'
            gtf = read.table("${gtf}", sep="\t")
            gtf = subset(gtf, V3 == "exon")
            write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
        RCODE

        awk '{print $1 "\t" ($2 - 1) "\t" $3}' exome.bed > exome.fixed.bed

        ${gatk_path} \
            BedToIntervalList \
            -I exome.fixed.bed \
            -O ${output_name} \
            -SD ${ref_dict}
    
    >>>

    output {
        File interval_list = "${output_name}"
    }

    runtime {
        docker: docker
        preemptible: preemptible_count
    }
}

#NOTE: assuming aggregated bams & paired end fastqs
task SamToFastq {
    input {
        File unmapped_bam
        String base_name

        String gatk_path

        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
            SamToFastq \
            --INPUT ${unmapped_bam} \
            --VALIDATION_STRINGENCY SILENT \
            --FASTQ ${base_name}.1.fastq.gz \
            --SECOND_END_FASTQ ${base_name}.2.fastq.gz
    }

    output {
        File fastq1 = "${base_name}.1.fastq.gz"
        File fastq2 = "${base_name}.2.fastq.gz"
    }

    runtime {
        docker: docker
        memory: "4 GB"
        disks: "local-disk " + ((size(unmapped_bam,"GB")+1)*5) + " HDD"
        preemptible: preemptible_count
    }
}

task StarGenerateReferences {
    input {
        File ref_fasta
        File ref_fasta_index
        File annotations_gtf
        Int? read_length  ## Should this be an input, or should this always be determined by reading the first line of a fastq input

        Int? num_threads
        Int threads = select_first([num_threads, 8])
        
        Int? additional_disk
            Int add_to_disk = select_first([additional_disk, 0])
        Int disk_size = select_first([100 + add_to_disk, 100])
        Int? mem_gb
        Int mem = select_first([100, mem_gb])
        String docker
        Int preemptible_count
    }

    command {
        set -e
        mkdir STAR2_5

        STAR \
        --runMode genomeGenerate \
        --genomeDir STAR2_5 \
        --genomeFastaFiles ${ref_fasta} \
        --sjdbGTFfile ${annotations_gtf} \
        ${"--sjdbOverhang "+(read_length-1)} \
        --runThreadN ${threads}

        ls STAR2_5

        tar -zcvf star-HUMAN-refs.tar.gz STAR2_5
    }

    output {
        Array[File] star_logs = glob("*.out")
        File star_genome_refs_zipped = "star-HUMAN-refs.tar.gz"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size + " HDD"
        cpu: threads
        memory: mem +" GB"
        preemptible: preemptible_count
    }
}


task StarAlign {
    input {
        File star_genome_refs_zipped
        File fastq1
        File fastq2
        String base_name
        Int? read_length

        Int? num_threads
        Int threads = select_first([num_threads, 8])
        Int? star_mem_max_gb
        Int star_mem = select_first([star_mem_max_gb, 45])
        #Is there an appropriate default for this?
        Int? star_limitOutSJcollapsed

        Int? additional_disk
        Int add_to_disk = select_first([additional_disk, 0])
        String docker
        Int preemptible_count
    }

    command <<<
        set -e

        tar -xvzf ${star_genome_refs_zipped}

        STAR \
        --genomeDir STAR2_5 \
        --runThreadN ${threads} \
        --readFilesIn ${fastq1} ${fastq2} \
        --readFilesCommand "gunzip -c" \
        ${"--sjdbOverhang "+(read_length-1)} \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --limitBAMsortRAM ${star_mem+"000000000"} \
        --limitOutSJcollapsed ${default=1000000 star_limitOutSJcollapsed} \
        --outFileNamePrefix ${base_name}.
    >>>

    output {
        File output_bam = "${base_name}.Aligned.sortedByCoord.out.bam"
        File output_log_final = "${base_name}.Log.final.out"
        File output_log = "${base_name}.Log.out"
        File output_log_progress = "${base_name}.Log.progress.out"
        File output_SJ = "${base_name}.SJ.out.tab"
    }

    runtime {
        docker: docker
        disks: "local-disk " + ((size(fastq1,"GB")+size(fastq2,"GB")*10)+30+add_to_disk) + " HDD"
        memory: (star_mem+1) + " GB"
        cpu: threads
        preemptible: preemptible_count
    }
}

task MergeBamAlignment {
    input {

        File ref_fasta
        File ref_dict

        File unaligned_bam
        File star_bam
        String base_name

        String gatk_path

        String docker
        Int preemptible_count
        #Using default for max_records_in_ram
    }

    command {
        ${gatk_path} \
            MergeBamAlignment \
            --REFERENCE_SEQUENCE ${ref_fasta} \
            --UNMAPPED_BAM ${unaligned_bam} \
            --ALIGNED_BAM ${star_bam} \
            --OUTPUT ${base_name}.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT
    }

    output {
        File output_bam="${base_name}.bam"
    }

    runtime {
        docker: docker
        disks: "local-disk " + ((size(unaligned_bam,"GB")+size(star_bam,"GB")+1)*5) + " HDD"
        memory: "4 GB"
        preemptible: preemptible_count
    }
}

task MarkDuplicates {
    input {

        File input_bam
        String base_name

        String gatk_path

        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
            MarkDuplicates \
            --INPUT ${input_bam} \
            --OUTPUT ${base_name}.bam  \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --METRICS_FILE ${base_name}.metrics
    }

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
        File metrics_file = "${base_name}.metrics"
    }

    runtime {
        disks: "local-disk " + ((size(input_bam,"GB")+1)*3) + " HDD"
        docker: docker
        memory: "4 GB"
        preemptible: preemptible_count
    }
}

task SplitNCigarReads {
    input {
        File input_bam
        File input_bam_index
        String base_name

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String gatk_path
        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
                SplitNCigarReads \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${base_name}.bam 
    }

        output {
                File output_bam = "${base_name}.bam"
                File output_bam_index = "${base_name}.bai"
        }

    runtime {
        disks: "local-disk " + ((size(input_bam,"GB")+1)*5 + size(ref_fasta,"GB")) + " HDD"
        docker: docker
        memory: "4 GB"
        preemptible: preemptible_count
    }
}

task BaseRecalibrator {
    input {

        File input_bam
        File input_bam_index
        String recal_output_file

        File dbSNP_vcf
        File dbSNP_vcf_index
        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String gatk_path

        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${recal_output_file} \
            -known-sites ${dbSNP_vcf} \
            -known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
    }

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
        memory: "6 GB"
        disks: "local-disk " + (size(input_bam,"GB")*3)+30 + " HDD"
        docker: docker
        preemptible: preemptible_count
    }
}


task ApplyBQSR {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File recalibration_report

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String gatk_path

        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${base_name}.bam \
            --bqsr-recal-file ${recalibration_report}
    }

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        memory: "3500 MB"
        disks: "local-disk " + (size(input_bam,"GB")*4)+30 + " HDD"
        preemptible: preemptible_count
        docker: docker
    }
}

task HaplotypeCaller {
    input {

        File input_bam
        File input_bam_index
        String base_name

        File interval_list

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File dbSNP_vcf
        File dbSNP_vcf_index

        String gatk_path
        String docker
        Int preemptible_count

        Int? stand_call_conf
    }

    command {
        ${gatk_path} --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -L ${interval_list} \
        -O ${base_name}.vcf.gz \
        -dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
        --dbsnp ${dbSNP_vcf}
    }

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        memory: "6.5 GB"
        disks: "local-disk " + (size(input_bam,"GB")*2)+30 + " HDD"
        preemptible: preemptible_count
    }
}

task VariantFiltration {
    input {

        File input_vcf
        File input_vcf_index
        String base_name

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String gatk_path
        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
            VariantFiltration \
            --R ${ref_fasta} \
            --V ${input_vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" \
            --filter "FS > 30.0" \
            --filter-name "QD" \
            --filter "QD < 2.0" \
            -O ${base_name}
    }

    output {
        File output_vcf = "${base_name}"
        File output_vcf_index = "${base_name}.tbi"
    }

    runtime {
        docker: docker
        memory: "3 GB"
        disks: "local-disk " + (size(input_vcf,"GB")*2)+30 + " HDD"
        preemptible: preemptible_count
    }
}

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String output_vcf_name

        Int disk_size = 5

        String gatk_path

        String docker
        Int preemptible_count

    # Using MergeVcfs instead of GatherVcfs so we can create indices
        }
        # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command {
        ${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT ' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    }

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    runtime {
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_count
    }
}

task ScatterIntervalList {
    input {

        File interval_list
        Int scatter_count
        String gatk_path
        String docker
        Int preemptible_count
    }

    command {
        set -e
        mkdir out
        ${gatk_path} --java-options "-Xms1g" \
            IntervalListTools \
            --SCATTER_COUNT ${scatter_count} \
            --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            --UNIQUE true \
            --SORT true \
            --INPUT ${interval_list} \
            --OUTPUT out
    
        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
            (directory, filename) = os.path.split(interval)
            newName = os.path.join(directory, str(i + 1) + filename)
            os.rename(interval, newName)
        print(len(intervals))
        if len(intervals) == 0:
            raise ValueError("Interval list produced 0 scattered interval lists. Is the gtf or input interval list empty?")
        f = open("interval_count.txt", "w+")
        f.write(str(len(intervals)))
        f.close()
        CODE
    }

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

    runtime {
        disks: "local-disk 1 HDD"
        memory: "2 GB"
        docker: docker
        preemptible: preemptible_count
    }
}

task RevertSam {
    input {
        File input_bam
        String base_name
        String sort_order

        String gatk_path

        String docker
        Int preemptible_count
    }

    command {
        ${gatk_path} \
            RevertSam \
            --INPUT ${input_bam} \
            --OUTPUT ${base_name}.bam \
            --VALIDATION_STRINGENCY SILENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --SORT_ORDER ${sort_order}
    }

    output {
        File output_bam = "${base_name}.bam"
    }

    runtime {
        docker: docker
        disks: "local-disk " + (size(input_bam,"GB")+1)*5 + " HDD"
        memory: "4 GB"
        preemptible: preemptible_count
    }
}

task Funcotate {
    input {
        File ref_fasta
        File ref_fai
        File ref_dict
        File input_vcf
        File input_vcf_idx
        String reference_version
        String output_file_base_name
        String output_format="VCF"
        Boolean compress=true
        Boolean use_gnomad=false
        # This should be updated when a new version of the data sources is released
        # TODO: Make this dynamically chosen in the command.
        File data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz"
        String? control_id
        String? case_id
        String? sequencing_center
        String? sequence_source
        String? transcript_selection_mode
        File? transcript_selection_list
        Array[String]? annotation_defaults
        Array[String]? annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean filter_funcotations=true
        File? interval_list

        String? extra_args
        String? gcs_project_for_requester_pays

        # ==============
        Int? disk_space   #override to request more disk than default small task params

        # You may have to change the following two parameter values depending on the task requirements
        Int default_ram_mb = 3000
        # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
        Int default_disk_space_gb = 100
        
        String gatk_docker
        Int boot_disk_size=12
        Int machine_mem=4000
        Int preemptible=2
        Int max_retries=2
        Int cpu=2
        }

        # ==============
        # Process input args:
        String output_maf = output_file_base_name + ".maf"
        String output_maf_index = output_maf + ".idx"
        String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
        String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
        String output_file = if output_format == "MAF" then output_maf else output_vcf
        String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
        String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
        String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
        String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
        String filter_funcotations_args = if filter_funcotations then " --remove-filtered-variants " else ""
        String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
        String interval_list_arg = if defined(interval_list) then " -L " else ""
        String extra_args_arg = select_first([extra_args, ""])

        String dollar = "$"

        parameter_meta{
            ref_fasta: {localization_optional: true}
            ref_fai: {localization_optional: true}
            ref_dict: {localization_optional: true}
            input_vcf: {localization_optional: true}
            input_vcf_idx: {localization_optional: true}
        }

        command <<<
            set -e
            export GATK_LOCAL_JAR="/root/gatk.jar"

            # Extract our data sources:
            echo "Extracting data sources zip file..."
            mkdir datasources_dir
            tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
            DATA_SOURCES_FOLDER="$PWD/datasources_dir"

            # Handle gnomAD:
            if ~{use_gnomad} ; then
                echo "Enabling gnomAD..."
                for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                    if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                        cd ~{dollar}{DATA_SOURCES_FOLDER}
                        tar -zvxf ~{dollar}{potential_gnomad_gz}
                        cd -
                    else
                        echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                        false
                    fi
                done
            fi

        # Run Funcotator:
        gatk --java-options "-Xmx~{machine_mem - 500}m" Funcotator \
            --data-sources-path $DATA_SOURCES_FOLDER \
            --ref-version ~{reference_version} \
            --output-file-format ~{output_format} \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            -O ~{output_file} \
            ~{interval_list_arg} ~{default="" interval_list} \
            --annotation-default normal_barcode:~{default="Unknown" control_id} \
            --annotation-default tumor_barcode:~{default="Unknown" case_id} \
            --annotation-default Center:~{default="Unknown" sequencing_center} \
            --annotation-default source:~{default="Unknown" sequence_source} \
            ~{"--transcript-selection-mode " + transcript_selection_mode} \
            ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
            ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
            ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
            ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
            ~{filter_funcotations_args} \
            ~{extra_args_arg} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
        # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
        if [[ "~{output_format}" == "MAF" ]] ; then
        touch ~{output_maf_index}
        fi
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File funcotated_output_file = "~{output_file}"
        File funcotated_output_file_index = "~{output_file_index}"
    }
}