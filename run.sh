### Build the Docker container
sudo docker build -t kauralasoo/coloc-wrapper .

### Push to DockerHub
docker push kauralasoo/coloc-wrapper

### Build a local copy of the Singularity container
singularity build coloc-wrapper.img docker://kauralasoo/coloc-wrapper:latest

### Bind the /gpfs/hpc file system to all singularity containers (helps with testing)
### You can add this line to your ~/.bashrc file
export SINGULARITY_BINDPATH="/gpfs/hpc"

### Run coloc all QTLs agains all GWAS traits
snakemake --use-singularity --cluster scripts/snakemake_submit_UT.py -np -s run_coloc.snakefile results/coloc/coloc_out.txt --jobs 300 --configfile config/coloc_config.yaml --rerun-incomplete