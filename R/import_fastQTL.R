
#' Fetch particular genes from tabix indexed FastQTL output file.
#'
#' @param phenotype_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing Rasqual results for each gene.
#' @export
fastqtlTabixFetchGenes <- function(phenotype_ranges, tabix_file){

  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))

  #Set column names for rasqual
  fastqtl_columns = c("phenotype_id","chr","pos","snp_id","distance","p_nominal","beta")
  fastqtl_coltypes = "ccicidd"

  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], col_names = fastqtl_columns, col_types = fastqtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)

    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}

#' Import fastQTL output table into R.
#'
#' Detect if the table is from nominal run or permutation run and add proper column names.
#'
#' @param file_path Path to the fastQTL output file
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export
importFastQTLTable <- function(file_path){
  table = read.table(file_path, stringsAsFactors = FALSE)
  if(ncol(table) == 11){
    colnames(table) = c("gene_id", "n_cis_snps", "beta1", "beta2", "dummy", "snp_id", "distance","p_nominal", "slope","p_perm","p_beta")
    table = table %>% tbl_df() %>%
      dplyr::filter(!is.na(p_beta)) %>%
      dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
      dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues) %>%
      dplyr::arrange(p_fdr)
  }
  return(table)
}
