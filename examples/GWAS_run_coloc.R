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
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("--qtl"), type="character", default=NULL,
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type"),
  make_option(c( "--gwasvarinfo"), type="character", default=NULL,
              help="Variant infromation file for the GWAS dataset.", metavar = "type"),
  make_option(c("--qtlvarinfo"), type="character", default=NULL,
              help="Variant information file for the QTL dataset.", metavar = "type"),
  make_option(c("--gwaslist"), type="character", default=NULL,
              help="Path to the list of GWAS studies.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
opt = list(gwas = "IBD", w = "2e5", p = "leafcutter", d = "~/datasets/Inflammatory_GWAS/", o = "results/acLDL/coloc/coloc_lists/", qtl = "processed/salmonella/qtltools/output/", s = "analysis/data/sample_lists/salmonella_coloc_sample_sizes.txt", gwasvarinfo = "results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz", qtlvarinfo = "results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz", gwaslist = "analysis/data/gwas/GWAS_summary_stat_list.labeled.txt")
opt$d = "/nfs/users/nfs_k/ka8/scratch/datasets/Inflammatory_GWAS/"

#Extract parameters for CMD options
gwas_id = opt$gwas
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
qtl_dir = opt$qtl
outdir = opt$o
sample_size_path = opt$s
gwas_var_path = opt$gwasvarinfo
qtl_var_path = opt$qtlvarinfo
gwas_list = opt$gwaslist

#Import variant information
gwas_var_info = colocWrapper::importVariantInformation(gwas_var_path)
qtl_var_info = colocWrapper::importVariantInformation(qtl_var_path)
print("Variant information imported.")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type"), col_type = "ccc")

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "ci")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name


#Construct a new QTL list
phenotype_values = colocWrapper::constructQtlListForColoc(phenotype, qtl_dir, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = colocWrapper::prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix,
                                       gwas_variant_info = gwas_var_info, fdr_thresh = 0.1,
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
print("Pre-filtering completed.")

#Test for coloc
coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes,
                             ~colocWrapper::colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x,
                                                      gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"),
                                                      gwas_variant_info = gwas_var_info,
                                                      qtl_variant_info = qtl_var_info,
                                                      N_qtl = .y, cis_dist = cis_window))
print("Coloc completed.")

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w, "txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
colocMolecularQTLs(qtl_pairs[1,], phenotype_values$qtl_summary_list$naive, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), gwas_var_info, qtl_var_info, N_qtl = 84, cis_dist = 2e5)
