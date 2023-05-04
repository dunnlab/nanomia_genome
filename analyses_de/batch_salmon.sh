#!/bin/bash

#SBATCH --job-name=salmon
#SBATCH --output=salmon.log
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16           # number of cores
#SBATCH --mem-per-cpu=5G             # shared memory, scaling with CPU request

# Set up paths
cwd=$(pwd)
echo $cwd
base=$(basename $cwd)
echo $base
scratch_dir="/gpfs/ysm/scratch60/dunn/cwd7/$base"
echo "Using this scratch dir for output: "
echo $scratch_dir
# mkdir -p $scratch_dir

# Set up modules
module purge # Unload any existing modules that might conflict
module load miniconda
module list

conda activate salmon

REF="PO2744_Nanomia_bijuga.transcript.fasta"
READ_DIR="/gpfs/gibbs/project/dunn/cwd7/nanomia_cat_reads/"

# from https://combine-lab.github.io/salmon/getting_started/

salmon index -t ${REF} -i ${REF}.index

mkdir -p quant

for fastq_file in ${READ_DIR}*.fastq; do

  filename=$(basename "$fastq_file")
  echo "Processing sample ${filename}"
  salmon quant -i ${REF}.index -l A --gcBias \
         -r ${fastq_file} \
         -p 16 -o quant/${filename}.quant
done 
