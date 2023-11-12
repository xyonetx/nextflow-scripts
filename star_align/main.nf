// This script simply creates the STAR index

process star_align {

    publishDir "${params.output_dir}/bams", mode:"copy", pattern: "*.bam"
    publishDir "${params.output_dir}/star_logs", mode:"copy", pattern: "*.out"
    publishDir "${params.output_dir}/quantifications", mode:"copy", pattern: "*.tsv"
    container "ghcr.io/xyonetx/nextflow-scripts/star:2.7.11a"
    cpus 8
    memory '32 GB'

    input:
        path star_index_tar
        tuple val(sampleID), path("reads_R?.fastq.gz")

    output:
        path("*.sorted.bam")
        path("*.quantifications.tsv")
        path("*.Log.final.out")

    script:
        """
        # Params to match those of the GDC mRNA analysis pipeline:
        # https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

        mkdir ./starIndex
        tar -xf ${star_index_tar} -C ./starIndex

        # First pass alignment
        /opt/software/STAR-2.7.11a/bin/Linux_x86_64_static/STAR \
            --runThreadN 8 \
            --twopassMode Basic \
            --genomeDir ./starIndex \
            --readFilesIn reads_R1.fastq.gz reads_R2.fastq.gz \
            --readFilesCommand zcat \
            --sjdbOverhang ${params.sjdb_overhang} \
            --outFilterMultimapScoreRange 1 \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 10 \
            --alignIntronMax 500000 \
            --alignMatesGapMax 1000000 \
            --sjdbScore 2 \
            --alignSJDBoverhangMin 1 \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterScoreMinOverLread 0.33 \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts

        mv Aligned.sortedByCoord.out.bam ${sampleID}.sorted.bam
        mv ReadsPerGene.out.tab ${sampleID}.quantifications.tsv
        mv Log.final.out ${sampleID}.Log.final.out
        """
}


process multiqc {

    cpus 2
    memory '4 GB'
    container "ghcr.io/xyonetx/nextflow-scripts/star:2.7.11a"
    publishDir "${params.output_dir}/multiqc", mode:"copy"

    input:
        path("*.Log.final.out")

    output:
        path('multiqc_report.html')
        path('multiqc_data', type: 'dir')

    script:
        """
        multiqc .
        """
}

workflow {    
    fq_channel = Channel.fromFilePairs(params.fastq_dir + '/' + params.fastq_pattern)
    (bams_ch, quants_ch, logs_ch) = star_align(params.star_index_path, fq_channel)
    multiqc(logs_ch.collect())
}