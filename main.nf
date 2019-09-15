Channel.fromPath(params.sumstats)
    .ifEmpty { error "Cannot find the sumstats file in: ${params.sumstats}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ file(row.sumstats_file)]}
    .set { sumstats_channel }

process index_gwas_catalog_sumstats{

    input:
    file(sumstats) from sumstats_channel

    output:
    set file("${sumstats.baseName}.tsv.gz"), file("${sumstats.baseName}.tsv.gz.tbi") to indexed_sumstats

    script:
    """
    zcat ${sumstats} | tail -n+2 | awk '{if($3 != "NA") print $0}' | LANG=C sort -k3,3 -k4,4n | bgzip > ${sumstats.baseName}.tsv.gz
    tabix -b4 -e4 -s3 ${sumstats.baseName}.tsv.gz
    """
}