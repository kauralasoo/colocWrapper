library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
devtools::load_all("colocWrapper/")

#Parse command-line options
option_list <- list(
  make_option(c("--qtl_leads"), type="character", default=NULL,
              help="Path to the QTL lead variants used for coloc.", metavar = "type"),
  make_option(c("--qtl_stats"), type="character", default=NULL,
              help="Path to the full QTL summary statistics file.", metavar = "type"),
  make_option(c("--cis_window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas_leads"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("--gwas_stats"), type="character", default=NULL,
              help="Path to GWAS summary stats file.", metavar = "type"),
  make_option(c("--out"), type="character", default=NULL,
              help="Path to the output file.", metavar = "type"),
  make_option(c("--qtl_varinfo"), type="character", default=NULL,
              help="Variant information file for the QTL dataset.", metavar = "type"),
  make_option(c("--gwas_type"), type="character", default=NULL,
              help="GWAS summary statistics type.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
if(FALSE){
  opt = list(qtl_leads = "testdata/platelet.permuted.txt.gz",
             qtl_stats = "testdata/platelet.nominal.sorted.txt.gz",
             cis_window = 200000,
             gwas_leads = "testdata/27863252-GCST004599-EFO_0004584.top_hits.tsv.gz",
             gwas_stats = "testdata/27863252-GCST004599-EFO_0004584.tsv.gz",
             qtl_varinfo = "testdata/platelet.variant_information.txt.gz",
             out = "colocalised_hits.txt",
             gwas_type = "GWASCatalog")
}


#Extract parameters for CMD options
cis_window = as.numeric(opt$cis_window)
outfile = opt$out
gwas_id = basename(opt$gwas_stats) %>%
  sub('\\.tsv.gz$', '', .)

#Import variant information
varinfo_df = colocWrapper::importVariantInformation(opt$qtl_varinfo)
print("Variant information imported.")

#Extract sample size from variant info
sample_size = median(varinfo_df$AN)/2

#Construct a new QTL list
qtl_list = colocWrapper::constructQtlmapListForColoc(lead_path = opt$qtl_leads, summary_path = opt$qtl_stats, sample_size)

#Prefilter coloc candidates
qtl_df_list = colocWrapper::prefilterColocCandidates(qtl_list$min_pvalues, opt$gwas_leads,
                                       variant_info = varinfo_df, fdr_thresh = 0.1,
                                       overlap_dist = 1e5, gwas_type = opt$gwas_type)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
message("Pre-filtering completed. Number of loci included for colocalisation: ", nrow(qtl_pairs))

#Test for coloc
coloc_res_list = purrr::map2(qtl_list$qtl_summary_list, qtl_list$sample_sizes,
                             ~colocWrapper::colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x,
                                                      gwas_summary_path = opt$gwas_stats,
                                                      variant_info = varinfo_df,
                                                      N_qtl = .y, cis_dist = cis_window, gwas_type = opt$gwas_type))
message("Coloc completed.")

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity) %>%
  dplyr::arrange(-PP.H4.abf) %>%
  dplyr::select(-.row) %>%
  dplyr::mutate(gwas_trait = gwas_id)
write.table(coloc_hits, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging
if(FALSE){
  colocMolecularQTLs(qtl_pairs, qtl_list$qtl_summary_list$qtl_group, gwas_summary_path = opt$gwas_stats,
                     variant_info = varinfo_df, N_qtl = qtl_list$sample_sizes$qtl_group, cis_dist = cis_window,
                     gwas_type = opt$gwas_type)
}
