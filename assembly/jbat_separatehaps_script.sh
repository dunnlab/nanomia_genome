#!/bin/bash
#SBATCH --job-name=jbat
#SBATCH --output=chromap.log
#SBATCH --time=8-00:00:00
#SBATCH --partition=pi_dunn
#SBATCH --nodes=1 --ntasks-per-node=32                    # number of cores and nodes
#SBATCH --mem=180GB                  # memory pool for all cores

# Max on pi_dunn
# Maximum job time limit	14-00:00:00
# 36 cores, 181 GB RAM

# Author: Namrata Ahuja
# Date  : August 2022
# email : namrata.ahuja@yale.edu 

PREFIX=NS_hap1_hap2_v0.8

THREAD=32


# Directory roots

PROJECT_DIR=/gpfs/ysm/project/dunn/na375/20220818_chromap
SCRATCH_DIR=/gpfs/ysm/scratch60/dunn/na375/20220818_chromap

# Assembly
ASSEM=/gpfs/ysm/project/dunn/na375/20220818_chromap/NS_hap1_hap2_v0.8.fasta




# Reads from /SAY/standard/cwd7-CC0522-FASEEB/data/cdunn/2022_dovetail_deliveries_nanomia/20220816_omni-C
OMNI_DIR=${SCRATCH_DIR}/omni-c
LANE1_R1=${OMNI_DIR}/DTG-OmniC-410-run1_R1_001.fastq.gz
LANE1_R2=${OMNI_DIR}/DTG-OmniC-410-run1_R2_001.fastq.gz
LANE2_R1=${OMNI_DIR}/DTG-OmniC-410-run2_R1_001.fastq.gz
LANE2_R2=${OMNI_DIR}/DTG-OmniC-410-run2_R2_001.fastq.gz

MYBIN=${PROJECT_DIR}/bin

# echo $MYBIN

${MYBIN}/chromap/chromap -i -r ${ASSEM} -o ${PREFIX}.index

awk -f ${MYBIN}/generate-assembly-file-from-fasta.awk ${ASSEM} > ${PREFIX}.assembly

module load miniconda
conda activate mapping
python ${MYBIN}/assembly-from-fasta.py -c ${ASSEM} py_${PREFIX}
conda deactivate
module unload miniconda
module load Python

for qual in 0
do
${MYBIN}/chromap/chromap -t ${THREAD} --preset hic -x ${PREFIX}.index -r ${ASSEM} -1 "${LANE1_R1},${LANE2_R1}" -2 "${LANE1_R2},${LANE2_R2}" -q ${qual} -o ${PREFIX}_${qual}.pairs
grep -v '#' ${PREFIX}_${qual}.pairs | \
        awk '{{if($6!="+") $6=16; else $6=0; if($7!="+") $7=16; else $7=0}} \
            $2<=$4{{print $6, $2, $3, 0, $7, $4, $5, 1, "1 - - 1  - - -" }} \
            $4<$2{{print $7, $4, $5, 0, $6, $2, $3, 1, "1 - - 1  - - -" }}' > ${PREFIX}_${qual}.longp
bash ${MYBIN}/3d-dna/visualize/run-assembly-visualizer.sh -p false ${PREFIX}.assembly ${PREFIX}_${qual}.longp || true
mv ${PREFIX}.hic ${PREFIX}_hic_${qual}.hic
done