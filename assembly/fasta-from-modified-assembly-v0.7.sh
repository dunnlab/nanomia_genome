#!/bin/bash
#SBATCH --job-name=fasta_make
#SBATCH --output=chromap.log
#SBATCH --time=8-00:00:00
#SBATCH --partition=pi_dunn
#SBATCH --nodes=1 --ntasks-per-node=32                    # number of cores and nodes
#SBATCH --mem=150GB                  # memory pool for all cores

# Max on pi_dunn
# Maximum job time limit	14-00:00:00
# 36 cores, 181 GB RAM

# Author: Namrata Ahuja
# Date  : August 2022
# email : namrata.ahuja@yale.edu

#this script takes a modified review.assembly file that is exported by juicebox and converts into a fasta file
#do this all in chromap directory
#need pysam- so activated conda env mapping

module load miniconda 
conda activate mapping

#first argument=artisanal file to use to make new fasta file from modified assembly, 
#second= -N to remove N's
#third=modified assembly to use
#fourth= name of original fasta file for og assembly
#fifth= prefix/name of new fasta file

/home/na375/project/20220818_chromap/bin/artisanal/bin/assembly-to-fasta -N py_NBI_hap1_hap2_v0.6.review.assembly NBI_hap1_hap2_v0.6.fasta NBI_hap1_hap2_unordered_v0.7 