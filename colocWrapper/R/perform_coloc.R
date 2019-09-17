#' Construct a GRanges obejct corresponding to a cis region around one variant.
#'
#' @param variant_df data frame with at least snp_id column
#' @param variant_information data.frame from importVariantInformation() function
#' @param cis_dist Number of basepairs upstream and downstream of the variant.
#'
#' @return GRanges
#' @export
constructVariantRanges <- function(variant_df, variant_information, cis_dist){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(variant_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "chr"))
  assertthat::assert_that(assertthat::has_name(variant_information, "pos"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% variant_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Add variant info to variant df
  var_df = dplyr::left_join(variant_df, var_info, by = "snp_id")

  #Make a ranges object
  var_ranges = var_df %>%
    dplyr::rename(seqnames = chr) %>%
    dplyr::mutate(start = pos - cis_dist, end = pos + cis_dist, strand = "*") %>%
    dataFrameToGRanges()

  return(var_ranges)
}

#' Test colocalisation between molecular QTL and GWAS summary stats
#'
#' @param qtl QTL summary stats (p_nominal, MAF, beta, snp_id)
#' @param gwas GWAS summary stats(beta, se, MAF, log_OR)
#' @param N_qtl Sample size of the QTL mapping study
#'
#' @return coloc.abf result object
#' @export
colocQtlGWAS <- function(qtl, gwas, N_qtl, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05){

  #Check that QTL df has all correct names
  assertthat::assert_that(assertthat::has_name(qtl, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl, "beta"))
  assertthat::assert_that(assertthat::has_name(qtl, "MAF"))
  assertthat::assert_that(assertthat::has_name(qtl, "p_nominal"))

  #Check that GWAS df has all correct names
  assertthat::assert_that(assertthat::has_name(gwas, "beta"))
  assertthat::assert_that(assertthat::has_name(gwas, "se"))
  assertthat::assert_that(assertthat::has_name(gwas, "snp_id"))
  assertthat::assert_that(assertthat::has_name(gwas, "log_OR"))
  assertthat::assert_that(assertthat::has_name(gwas, "MAF"))

  #Count NAs for log_OR and beta
  log_OR_NA_count = length(which(is.na(gwas$log_OR)))
  beta_NA_count = length(which(is.na(gwas$beta)))

  #Remove GWAS SNPs with NA std error
  gwas = dplyr::filter(gwas, !is.na(se))

  #If beta is not specified then use log_OR
  if(beta_NA_count <= log_OR_NA_count){
    coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal,
                                                 N = N_qtl,
                                                 MAF = qtl$MAF,
                                                 type = "quant",
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id),
                                 dataset2 = list(beta = gwas$beta,
                                                 varbeta = gwas$se^2,
                                                 type = "cc",
                                                 snp = gwas$snp_id,
                                                 s = 0.5, #This is acutally not used, because we already specified varbeta above.
                                                 MAF = gwas$MAF),
                                 p1 = p1, p2 = p2, p12 = p12)
  } else{
    coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal,
                                                 N = N_qtl,
                                                 MAF = qtl$MAF,
                                                 type = "cc",
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id),
                                 dataset2 = list(beta = gwas$log_OR,
                                                 varbeta = gwas$se^2,
                                                 type = "cc",
                                                 snp = gwas$snp_id,
                                                 s = 0.5, #This is acutally not used, because we already specified varbeta above.
                                                 MAF = gwas$MAF),
                                 p1 = p1, p2 = p2, p12 = p12)
  }

  return(coloc_res)
}


summaryReplaceCoordinates <- function(summary_df, variant_information){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% summary_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }

  #Add new coordinates
  new_coords = dplyr::select(summary_df, -chr, -pos) %>%
    dplyr::left_join(var_info, by = "snp_id") %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::arrange(pos)

  return(new_coords)
}

summaryReplaceSnpId <- function(summary_df, variant_information){

  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))

  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, pos %in% summary_df$pos) %>%
    dplyr::select(snp_id, chr, pos, MAF)

  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }

  #Add new coordinates
  new_coords = dplyr::select(summary_df, -snp_id) %>%
    dplyr::left_join(var_info, by = c("chr","pos")) %>%
    dplyr::filter(!is.na(snp_id)) %>%
    dplyr::arrange(pos)

  return(new_coords)
}

#' Perform colocalisation between an single QTL and GWAS summary stats from the same region.
#'
#' @param qtl_df A data frame with a single row and two columns (phenotype_id, snp_id) corresponding the a single QTL.
#' @param qtl_summary_path Path to the tabix indexed QTL summary statistics file generated by QTLtools
#' @param gwas_summary_path Path to the tabix indexed GWAS summary statistics file.
#' @param variant_info Variant information for the QTL data. Imported using importVariantInformation() function.
#' @param N_qtl Sample size for QTL mapping. Used by coloc to estimate the the standard errors from p-values and effect sizes.
#' @param cis_dist With of the genomic region around the lead QTL variant that is used for colocalisation; width = 2*cis_dist.
#' @param gwas_type type of GWAS summary stats used for coloc.
#'
#' @return List of colocalisation results or NULL values if there was an error.
#' @export
colocMolecularQTLs <- function(qtl_df, qtl_summary_path, gwas_summary_path,
                               variant_info,
                               N_qtl = 84, cis_dist = 1e5, gwas_type){

  #Assertions
  assertthat::assert_that(assertthat::has_name(qtl_df, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "snp_id"))
  assertthat::assert_that(nrow(qtl_df) == 1)

  assertthat::assert_that(is.numeric(cis_dist))
  assertthat::assert_that(is.numeric(N_qtl))

  #Print for debugging
  print(qtl_df$phenotype_id)

  result = tryCatch({
    #Make GRanges object to fetch data
    qtl_ranges = constructVariantRanges(qtl_df, variant_info, cis_dist = cis_dist)

    #Fetch QTL summary stats
    qtl_summaries = qtltoolsTabixFetchPhenotypes(qtl_ranges, qtl_summary_path)[[1]] %>%
      dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta)

    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWAS(qtl_ranges, gwas_summary_path, gwas_type)[[1]]

    #Substitute coordinate for the eqtl summary stats and add MAF
    qtl = summaryReplaceCoordinates(qtl_summaries, variant_info)

    #Substitute snp_id for the GWAS summary stats and add MAF
    gwas = summaryReplaceSnpId(gwas_summaries, variant_info)

    #Extract minimal p-values for both traits
    qtl_min = dplyr::arrange(qtl, p_nominal) %>% dplyr::filter(row_number() == 1)
    gwas_min = dplyr::arrange(gwas, p_nominal) %>% dplyr::filter(row_number() == 1)

    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl, gwas, N_qtl = N_qtl)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_min$p_nominal, gwas_pval = gwas_min$p_nominal,
                    qtl_lead = qtl_min$snp_id, gwas_lead = gwas_min$snp_id) #Add minimal pvalues

    #Summary list
    data_list = list(qtl = qtl, gwas = gwas)

    result = list(summary = coloc_summary, data = data_list)
  }, error = function(err) {
    print(paste("ERROR:",err))
    result = list(summary = NULL, data = NULL)
  }
  )
  return(result)
}

#' Applies the colocMolecularQTLs function to each row of the qtl_df data frame.
#'
#' See documentation for colocMolecularQTLs for more details
#' @param qtl_df Data frame of QTLs
#' @param ... Additional parameters passed on to colocMolecularQTLs function
#'
#' @export
colocMolecularQTLsByRow <- function(qtl_df, qtl_summary_path, gwas_summary_path,
                                    variant_info,
                                    N_qtl, cis_dist, gwas_type){
  result = purrrlyr::by_row(qtl_df, ~colocMolecularQTLs(., qtl_summary_path, gwas_summary_path,
                                                        variant_info,
                                                        N_qtl, cis_dist, gwas_type)$summary, .collate = "rows")
}

#' Perform a quick pre-filtering between QTLs and GWAS hits to reduce the number of coloc tests
#'
#' @param qtl_min_pvalues List of data frames with QTL lead pvalues. Each data frame must contain
#' gene_id, snp_id and p_fdr and should not contain other columns.
#' @param gwas_leads_path Prefix of the GWAS summarystats file
#' @param variant_info QTL variant information in GWAS coordinates.
#' @param fdr_thresh Minimal QTL FDR threshold
#' @param overlap_dist Max distance between GWAS and QTL variants.
#'
#' @return List of data.frames with phenotype_ids and snp_ids to be tested with coloc.
#' @export
prefilterColocCandidates <- function(qtl_min_pvalues, gwas_leads_path, variant_info,
                                     fdr_thresh = 0.1, overlap_dist = 1e5, gwas_type){

  #Make sure that the qtl_df has neccessary columns
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "p_fdr"))
  assertthat::assert_that(ncol(qtl_min_pvalues[[1]]) == 3)


  #Import top GWAS p-values
  gwas_pvals = importGWAS(gwas_leads_path, gwas_type) %>%
    dplyr::transmute(chr = chr, gwas_pos = pos)

  #Filter lead variants
  qtl_hits = purrr::map(qtl_min_pvalues, ~dplyr::filter(., p_fdr < fdr_thresh))
  lead_variants = purrr::map_df(qtl_hits, identity) %>% unique()
  selected_variants = dplyr::filter(variant_info, snp_id %in% lead_variants$snp_id) %>%
    dplyr::select(chr, pos, snp_id)

  #Add GRCh37 coordinates
  qtl_pos = purrr::map(qtl_hits, ~dplyr::left_join(., selected_variants, by = "snp_id") %>%
                         dplyr::filter(!is.na(pos)))

  #Identify genes that have associated variants nearby (ignoring LD)
  qtl_df_list = purrr::map(qtl_pos, ~dplyr::left_join(., gwas_pvals, by = "chr") %>%
                             dplyr::mutate(distance = abs(gwas_pos - pos)) %>%
                             dplyr::filter(distance < overlap_dist) %>%
                             dplyr::select(phenotype_id, snp_id) %>% unique())

}

constructQtlmapListForColoc <- function(lead_path, summary_path, sample_size){
  #Import min pvalues
  min_pvalues = list()
  min_pvalues[["qtl_group"]] = importQTLtoolsTable(lead_path) %>% dplyr::select(phenotype_id, snp_id, p_fdr)

  #Set summary path
  qtl_summary_list = list()
  qtl_summary_list[["qtl_group"]] = summary_path

  #Sample size list
  sample_size_list = list()
  sample_size_list[["qtl_group"]] = sample_size

  return(list(min_pvalues = min_pvalues, qtl_summary_list = qtl_summary_list, sample_sizes = sample_size_list))
}


#Utility functions that are often shared between multiple pacakages, but are too small to deserve their own package.

#' Convert a data frame into a GRanges object
#'
#' Seqnames, strand, start and end columns are used as corresponding elements
#' in the GRanges object. Remaining columns are added into the elementMetadata data frame.
#'
#' @param df Input data frame (required columns: seqnames, start, end, strand)
#'
#' @return GRanges object construct from df.
#' @export
dataFrameToGRanges <- function(df){
  #Convert a data.frame into a GRanges object

  gr = GenomicRanges::GRanges(seqnames = df$seqnames,
                              ranges = IRanges::IRanges(start = df$start, end = df$end),
                              strand = df$strand)

  #Add metadata
  meta = dplyr::select(df, -start, -end, -strand, -seqnames)
  GenomicRanges::elementMetadata(gr) = meta

  return(gr)
}
