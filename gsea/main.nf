process prep_gct_and_cls_files {

    tag "Prep GCT and CLS files for GSEA"
    publishDir "${params.output_dir}/gsea", mode:"copy"
    container "ghcr.io/xyonetx/nextflow-scripts/gsea:4.3.2"
    cpus 2
    memory '4 GB'

    input:
        path(norm_counts)
        path(annotations)

    output:
        path("${gct_file}")
        path("${cls_file}")

    script:
        gct_file = "input.gct"
        cls_file = "input.cls"
        """
        /usr/bin/python3 /opt/software/scripts/prep_files.py \
            -f ${norm_counts} \
            -a ${annotations} \
            -g ${gct_file} \
            -c ${cls_file} \
            -t ${params.min_reads}
       """
}


process prep_rnk_files {

    tag "Prep RNK files for GSEA"
    publishDir "${params.output_dir}/gsea", mode:"copy"
    container "ghcr.io/xyonetx/nextflow-scripts/gsea:4.3.2"
    cpus 2
    memory '4 GB'

    input:
        path(dge_results)

    output:
        path("${rnk_file}")

    script:
        rnk_file = "${dge_results.baseName}.rnk"
        """    
        /usr/bin/python3 /opt/software/scripts/create_rnk_file.py \
            -f ${dge_results} \
            -o ${rnk_file}
       """
}


process run_gsea {

    tag "Run GSEA"
    publishDir "${params.output_dir}/gsea", mode:"copy"
    container "ghcr.io/xyonetx/nextflow-scripts/gsea:4.3.2"
    cpus 2
    memory '4 GB'

    input:
        path(gct_file)
        path(cls_file)
        tuple val(base_condition), val(experimental_condition)

    output:
        path("${contrast}.gsea_results.zip")

    script:
        def contrast_template = "%s_versus_%s"
        contrast = String.format(contrast_template, experimental_condition, base_condition)
        """

        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEA \
            -res "${gct_file}" \
            -cls "${cls_file}#${contrast}" \
            -gmx  /opt/software/resources/h.all.v2023.2.Hs.symbols.gmt \
            -chip /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -out /gsea/ \
            -rpt_label output \
            -zip_report true \
            -collapse Collapse \
            -mode Max_probe \
            -norm meandiv \
            -nperm 1000 \
            -permute phenotype \
            -rnd_seed timestamp \
            -rnd_type no_balance \
            -scoring_scheme weighted \
            -metric Signal2Noise \
            -sort real \
            -order descending \
            -create_gcts false \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -median false \
            -num 100 \
            -plot_top_x 20 \
            -save_rnd_lists false \
            -set_max 500 \
            -set_min 15

        /usr/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea/output*/*.zip" \
            -o ${contrast}.gsea_results.zip
        """
}


process run_gsea_preranked {

    tag "Run GSEA Preranked"
    publishDir "${params.output_dir}/gsea_preranked", mode:"copy"
    container "ghcr.io/xyonetx/nextflow-scripts/gsea:4.3.2"
    cpus 2
    memory '4 GB'

    input:
        path(rnk_file)
        tuple val(base_condition), val(experimental_condition)

    output:
        path("${contrast}.gsea_preranked_results.zip")

    script:
        def contrast_template = "%s_versus_%s"
        contrast = String.format(contrast_template, experimental_condition, base_condition)
        """
        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
            -gmx /opt/software/resources/h.all.v2023.2.Hs.symbols.gmt \
            -chip /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -rnk ${rnk_file} \
            -rpt_label output \
            -collapse Collapse \
            -mode Abs_max_of_probes \
            -norm meandiv \
            -nperm 1000 \
            -rnd_seed timestamp \
            -scoring_scheme weighted \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -plot_top_x 20 \
            -set_max 500 \
            -set_min 15 \
            -zip_report true \
            -out /gsea_preranked/

        /usr/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea_preranked/output*/*.zip" \
            -o ${contrast}.gsea_preranked_results.zip
        """
}

workflow {
    dge_results_ch = Channel.fromPath(params.dge_results_glob)
    norm_counts_ch= Channel.fromPath(params.norm_counts)
    ann_ch = Channel.fromPath(params.annotations)

    (gct_ch, cls_ch) = prep_gct_and_cls_files(norm_counts_ch, ann_ch)
    gct_ch = gct_ch.collect()
    cls_ch = cls_ch.collect()

    rnk_ch = prep_rnk_files(dge_results_ch)

    contrast_ch = Channel.fromPath(params.contrasts)
           .splitCsv(header: ['base_condition', 'experimental_condition'], skip: 1 )
           .map{
               row -> tuple(row.base_condition, row.experimental_condition)
           }
    
    run_gsea(gct_ch, cls_ch, contrast_ch)
    run_gsea_preranked(rnk_ch, contrast_ch)
}