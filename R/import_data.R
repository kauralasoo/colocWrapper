#' Import variant information extracted from VCF file into R
#'
#' The variant information text file can be generated from the VCF using the following
#' bcftools command:
#' bcftools query -f '\%CHROM\\t\%POS\\t\%ID\\t\%REF\\t\%ALT\\t\%TYPE\\t\%AC\\t\%AN\\n' path/to/vcf_file.vcf.gz | bgzip > path/to/variant_infromation_file.txt.gz
#'
#' @param path Path to the the variant information text file.
#' @export
importVariantInformation <- function(path){
  info_col_names = c("chr","pos","snp_id","ref","alt","type","AC","AN")
  into_col_types = "ciccccii"
  snp_info = readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
  snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
    dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE)) %>%
    dplyr::mutate(MAF = pmin(AC/AN, 1-(AC/AN)))
  return(snp_info)
}

#' Import a specific region from a tabix-indexed GWAS summary stats file
tabixFetchGWASSummary <- function(granges, summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}

#' Import full GWAS summary stats file
importGWASSummary <- function(summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvals = readr::read_tsv(summary_path,
                               col_names = gwas_col_names, col_types = gwas_col_types, skip = 1)
  return(gwas_pvals)
}


#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}


