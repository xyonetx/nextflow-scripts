process deseq2_dge {

    cpus 4
    memory '12 GB'
    container "ghcr.io/xyonetx/nextflow-scripts/deseq2:1.40.2"
    publishDir "${params.output_dir}/diffential_expression/${contrast}", mode:"copy"

    input:
        path ann
        path raw_counts
        tuple val(base_condition), val(experimental_condition)

    output:
        path("${output_dge_filename}")
        path("normalized_counts.tsv")

    script:
        def output_dge_filename_template = "deseq2_results.%s.tsv"
        def contrast_template = "%s_vs_%s"
        contrast = String.format(contrast_template, experimental_condition, base_condition)
        output_dge_filename = String.format(output_dge_filename_template, contrast)
        """
        /opt/software/deseq2_env/bin/Rscript /opt/software/scripts/deseq2.R \
            ${raw_counts} \
            ${ann} \
            ${base_condition} \
            ${experimental_condition} \
            ${output_dge_filename} \
            normalized_counts.tsv
        """
}


process map_ensg_to_symbol {
    tag "Run ENSG to symbol gene mapping"
    publishDir "${output_dir}/${output_location}", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 2
    memory '4 GB'

    input:
        path infile
        val output_location

    output:
        path("${outfile}")

    script:
        outfile = "${infile.baseName}.symbol_remapped.tsv"
        """
        /usr/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${infile} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o ${outfile}
        """
}


workflow {    
    counts_ch= Channel.fromPath(params.raw_counts).collect()
    ann_ch = Channel.fromPath(params.annotations).collect()
    contrast_ch = Channel.fromPath(params.contrasts)
           .splitCsv(header: ['base_condition', 'experimental_condition'], skip: 1 )
           .map{
               row -> tuple(row.base_condition, row.experimental_condition)
           }

    (dge_ch, nc_ch) = deseq2_dge(ann_ch, counts_ch, contrast_ch)
    map_ensg_to_symbol(dge_ch, "differential_expression")
    map_ensg_to_symbol(nc_ch, "differential_expression")
}