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


#' Import full GWAS summary stats file
importGWASSummary <- function(summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvals = readr::read_tsv(summary_path,
                               col_names = gwas_col_names, col_types = gwas_col_types, skip = 1)
  return(gwas_pvals)
}

#' Import GWAS Catalog summary stats
importGWASCatalogSummary <- function(summary_path){
  gwas_col_names = c("hm_variant_id","hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta",
                "hm_odds_ratio","hm_ci_lower","hm_ci_upper","hm_effect_allele_frequency",
                "hm_code","variant","variant_id","chromosome","base_pair_location",
                "other_allele","effect_allele","alt_minor","direction","beta",
                "standard_error","p_value","mlog10p","effect_allele_frequency",
                "ma_freq","ci_lower","ci_upper","odds_ratio")
  gwas_pvals = readr::read_tsv(summary_path, col_names = gwas_col_names) %>%
    dplyr::transmute(snp_id = hm_variant_id, chr = as.character(hm_chrom), pos = hm_pos,
                     beta = hm_beta, OR = hm_odds_ratio, se = standard_error,
                     effect_AF = hm_effect_allele_frequency, rsid = hm_rsid, p_nominal = p_value) %>%
    dplyr::mutate(MAF = pmin(effect_AF, 1-effect_AF), log_OR = log(OR))
  return(gwas_pvals)
}

importGWAS <- function(summary_path, gwas_type = "GWASCatalog"){
  assertthat::assert_that(gwas_type %in% c("GWASCatalog", "Alasoo_2018"))

  if(gwas_type == "GWASCatalog"){
    result = importGWASCatalogSummary(summary_path)
    return(result)
  } else if (gwas_type == "Alasoo_2018"){
    result = importGWASSummary(summary_path)
    return(result)
  }
}

#' Import a specific region from a tabix-indexed GWAS summary stats file
tabixFetchGWASCatalogSummary <- function(granges, summary_path){
  gwas_col_names = c("hm_variant_id","hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta",
                     "hm_odds_ratio","hm_ci_lower","hm_ci_upper","hm_effect_allele_frequency",
                     "hm_code","variant","variant_id","chromosome","base_pair_location",
                     "other_allele","effect_allele","alt_minor","direction","beta",
                     "standard_error","p_value","mlog10p","effect_allele_frequency",
                     "ma_freq","ci_lower","ci_upper","odds_ratio")
  gwas_pvalues_list = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names)
  gwas_pvalues = purrr::map(gwas_pvalues_list, ~dplyr::transmute(.,snp_id = hm_variant_id,
                     chr = as.character(hm_chrom), pos = hm_pos,
                     beta = hm_beta, OR = hm_odds_ratio, se = standard_error,
                     effect_AF = hm_effect_allele_frequency, rsid = hm_rsid, p_nominal = p_value) %>%
    dplyr::mutate(MAF = pmin(effect_AF, 1-effect_AF), log_OR = log(OR)))
  return(gwas_pvalues)
}

#' Import a specific region from a tabix-indexed GWAS summary stats file
tabixFetchGWASSummary <- function(granges, summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF",
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}

tabixFetchGWAS <- function(granges, summary_path, gwas_type = "GWASCatalog"){
  assertthat::assert_that(gwas_type %in% c("GWASCatalog", "Alasoo_2018"))

  if(gwas_type == "GWASCatalog"){
    result = tabixFetchGWASCatalogSummary(granges, summary_path)
    return(result)
  } else if (gwas_type == "Alasoo_2018"){
    result = tabixFetchGWASSummary(granges, summary_path)
    return(result)
  }
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


