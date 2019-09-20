# colocWrapper
Simple workflow to perform colocalisation analysis between GWAS summary statistics from the GWASCatalog and eQTL summary statistics from the eQTL Catalogue.

## Example
For GWAS Catalog summary statistics, please use:
```bash
nextflow run main.nf\
 --gwasFile testdata/gwas_file.tsv\
 --studyFile testdata/study_file.tsv\
 --gwas_type 'GWASCatalog'\
 -profile coloc\
```

For internal summary statistics, please use:
```bash
nextflow run main.nf\
 --gwasFile testdata/gwas_file_Alasoo_2018.tsv\
 --studyFile testdata/study_file.tsv\
 --gwas_type 'Alasoo_2018'\
 -profile coloc
```

## Depedencies
The Docker container with all of the required dependencies to run the colocalisation pipeline is available from here:
https://hub.docker.com/r/kauralasoo/coloc-wrapper
