//sumstats_file
Channel.fromPath(params.gwasFile)
    .ifEmpty { error "Cannot find the sumstats file in: ${params.gwasFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ file(row.sumstats_file)]}
    .set { sumstats_channel }

// study	qtl_group	expression_matrix	phenotype_meta	sample_meta	vcf	phenotype_list	covariates
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find any qtl_results file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.study, row.qtl_group, row.quant_method, file(row.qtl_leads), file(row.qtl_stats), file("${row.qtl_stats}.tbi"), file(row.qtl_varinfo)]}
    .set { qtl_results_ch }

process index_gwas_catalog_sumstats{

    input:
    file(sumstats) from sumstats_channel

    output:
    set file("${sumstats.simpleName}.tsv.gz"), file("${sumstats.simpleName}.tsv.gz.tbi"), file("${sumstats.simpleName}.top_hits.tsv.gz") into indexed_sumstats

    script:
    """
    zcat ${sumstats} | tail -n+2 | awk '{if(\$3 != "NA") print \$0}' | LANG=C sort -k3,3 -k4,4n | bgzip > ${sumstats.simpleName}.tsv.gz
    tabix -b4 -e4 -s3 ${sumstats.simpleName}.tsv.gz
    zcat ${sumstats.simpleName}.tsv.gz | awk '{if(\$24 > ${params.gwas_min_logp}) print \$0}' | gzip > ${sumstats.simpleName}.top_hits.tsv.gz
    """
}

process run_coloc{

    input:
    set study, qtl_group, quant_method, file(qtl_leads), file(qtl_stats), file(qtl_stats_index), file(qtl_varinfo), file(gwas_stats), file(gwas_stats_index), file(gwas_leads) from qtl_results_ch.combine(indexed_sumstats)

    output:
    file("${study}.${qtl_group}.${quant_method}.${gwas_stats.baseName}.txt") into coloc_results

    script:
    """
    Rscript $baseDir/bin/qtlmap_run_coloc.R\
     --qtl_leads ${qtl_leads}\
     --qtl_stats ${qtl_stats}\
     --cis_window ${params.cis_window}\
     --gwas_leads ${gwas_leads}\
     --gwas_stats ${gwas_stats}\
     --out "${study}.${qtl_group}.${quant_method}.${gwas_stats.baseName}.txt"\
     --qtl_varinfo ${qtl_varinfo}\
     --gwas_type ${params.gwas_type}\
     --pkg_path $baseDir/colocWrapper
    """
}