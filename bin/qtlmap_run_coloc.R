suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("optparse"))

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
              help="GWAS summary statistics type.", metavar = "type"),
  make_option(c("--pkg_path"), type="character", default="colocWrapper/",
              help="GWAS summary statistics type.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
if(FALSE){
  opt = list(qtl_leads = "test_data/platelet.permuted.txt.gz",
             qtl_stats = "test_data/platelet.nominal.sorted.txt.gz",
             cis_window = 200000,
             gwas_leads = "test_data/Inflammatory_bowel_disease_UC_Liu_2015_NatGen_Immunochip.top_hits.GRCh38.txt.gz",
             gwas_stats = "test_data/Inflammatory_bowel_disease_UC_Liu_2015_NatGen_Immunochip.sorted.GRCh38.txt.gz",
             qtl_varinfo = "test_data/platelet.variant_information.txt.gz",
             out = "colocalised_hits.txt",
             gwas_type = "Alasoo_2018",
             pkg_path = "colocWrapper/")
}

#Load colocWrapper
devtools::load_all(opt$pkg_path)

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
  colocMolecularQTLs(qtl_pairs[1,], qtl_list$qtl_summary_list$qtl_group, gwas_summary_path = opt$gwas_stats,
                     variant_info = varinfo_df, N_qtl = qtl_list$sample_sizes$qtl_group, cis_dist = cis_window,
                     gwas_type = opt$gwas_type)
}
