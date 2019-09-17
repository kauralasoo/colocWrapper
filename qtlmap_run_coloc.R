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
opt = list(qtl_leads = "testdata/platelet.permuted.txt.gz",
            qtl_stats = "testdata/platelet.nominal.sorted.txt.gz",
            cis_window = 200000,
            gwas_leads = "testdata/27863252-GCST004599-EFO_0004584.top_hits.tsv.gz",
            gwas_stats = "testdata/27863252-GCST004599-EFO_0004584.tsv.gz",
            qtl_varinfo = "testdata/platelet.variant_information.txt.gz",
            out = "colocalised_hits.txt",
            gwas_type = "GWASCatalog")

#Extract parameters for CMD options
cis_window = as.numeric(opt$cis_window)
lead_vars = opt$qtl_leads
summary_path = opt$qtl_stats
outfile = opt$out
qtl_var_path = opt$qtl_varinfo

#Import variant information
qtl_var_info = colocWrapper::importVariantInformation(qtl_var_path)
print("Variant information imported.")

#Extract sample size from variant info
sample_size = median(qtl_var_info$AN)/2

#Import list of GWAS studies
#gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type"), col_type = "ccc")

#Construct a new QTL list
qtl_list = colocWrapper::constructQtlmapListForColoc(lead_path = lead_vars, summary_path = summary_path, sample_size)

#Spcecify the location of the GWAS summary stats file
#gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
#gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = colocWrapper::prefilterColocCandidates(qtl_list$min_pvalues, opt$gwas_leads,
                                       variant_info = qtl_var_info, fdr_thresh = 0.1,
                                       overlap_dist = 1e5, gwas_type = opt$gwas_type)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
print("Pre-filtering completed.")

#Test for coloc
coloc_res_list = purrr::map2(qtl_list$qtl_summary_list, qtl_list$sample_sizes,
                             ~colocWrapper::colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x,
                                                      gwas_summary_path = opt$gwas_stats,
                                                      variant_info = qtl_var_info,
                                                      N_qtl = .y, cis_dist = cis_window, gwas_type = opt$gwas_type))
print("Coloc completed.")

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity) %>%
  dplyr::arrange(-PP.H4.abf) %>%
  dplyr::mutate(gwas_trait = gwas_id)
write.table(coloc_hits, outfile, sep = "\t", quote = FALSE, row.names = FALSE)


#Debugging
if(FALSE){
  colocMolecularQTLs(qtl_pairs, qtl_list$qtl_summary_list$qtl_group, gwas_summary_path = opt$gwas_stats,
                     variant_info = qtl_var_info, N_qtl = qtl_list$sample_sizes$qtl_group, cis_dist = cis_window,
                     gwas_type = opt$gwas_type)
}
