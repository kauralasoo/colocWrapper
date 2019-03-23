### Build the Docker container
sudo docker build -t kauralasoo/coloc-wrapper .

### Push to DockerHub
docker push kauralasoo/coloc-wrapper

### Run coloc all QTLs agains all GWAS traits
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_coloc.snakefile results/coloc/coloc_out.txt --jobs 300 --configfile configs/coloc_config.yaml --rerun-incomplete