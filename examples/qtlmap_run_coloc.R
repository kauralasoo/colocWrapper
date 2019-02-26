library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../colocWrapper")

#Parse command-line options
option_list <- list(
  make_option(c("-l", "--leads"), type="character", default=NULL,
              help="File containing the lead variants used for colocalisation testing.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("--qtl"), type="character", default=NULL,
              help="Path to the full QTL summary statistics file.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("--qtlvarinfo"), type="character", default=NULL,
              help="Variant information file for the QTL dataset.", metavar = "type"),
  make_option(c("--gwaslist"), type="character", default=NULL,
              help="Path to the list of GWAS studies.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
#opt = list(gwas = "RA", w = "2e5",
#           l = "~/projects/RNAseq_pipeline/results/qtl_summary_stats/Fairfax_2014/array/monocyte_LPS2.permuted.txt.gz",
#           d = "~/datasets/GWAS_GRCh38/",
#           o = "results/monocyte_LPS2.coloc_results.txt",
#           qtl = "~/projects/RNAseq_pipeline/results/qtl_summary_stats/Fairfax_2014/array/monocyte_LPS2.nominal.sorted.txt.gz",
#           qtlvarinfo = "~/projects/RNAseq_pipeline/results/qtl_summary_stats/Fairfax_2014/array/monocyte_LPS2.variant_information.txt.gz",
#           gwaslist = "~/projects/macrophage-trQTLs/analysis/data/gwas/GWAS_summary_stat_list.labeled.txt")

#Extract parameters for CMD options
gwas_id = opt$gwas
cis_window = as.numeric(opt$w)
lead_vars = opt$l
gwas_dir = opt$d
summary_path = opt$qtl
outfile = opt$o
qtl_var_path = opt$qtlvarinfo
gwas_list = opt$gwaslist

#Import variant information
qtl_var_info = colocWrapper::importVariantInformation(qtl_var_path)
print("Variant information imported.")

#Extract sample size from variant info
sample_size = median(qtl_var_info$AN)/2

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type"), col_type = "ccc")

#Construct a new QTL list
qtl_list = colocWrapper::constructQtlmapListForColoc(lead_path = lead_vars, summary_path = summary_path, sample_size)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = colocWrapper::prefilterColocCandidates(qtl_list$min_pvalues, gwas_prefix,
                                       gwas_variant_info = qtl_var_info, fdr_thresh = 0.1,
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
print("Pre-filtering completed.")

#Test for coloc
coloc_res_list = purrr::map2(qtl_list$qtl_summary_list, qtl_list$sample_sizes,
                             ~colocWrapper::colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x,
                                                      gwas_summary_path = paste0(gwas_prefix, ".sorted.GRCh38.txt.gz"),
                                                      gwas_variant_info = qtl_var_info,
                                                      qtl_variant_info = qtl_var_info,
                                                      N_qtl = .y, cis_dist = cis_window))
print("Coloc completed.")

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity) %>%
  dplyr::arrange(-PP.H4.abf) %>%
  dplyr::mutate(gwas_trait = gwas_id)
write.table(coloc_hits, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], phenotype_values$qtl_summary_list$naive, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), qtl_var_info, qtl_var_info, N_qtl = 84, cis_dist = 2e5)
