#!/bin/bash
#
#SBATCH --job-name=octopus_v2
#SBATCH --cpus-per-task=20

# Author: Dalila Destanovic
# Date  : August 2022
# email : a12110418@unet.univie.ac.at or dalila.destanovic@gmail.com

PREFIX=octopus_vulgaris_v3
THREAD=20
ASSEM=octopus_vulgaris_v2_new_assembly.fasta
R1=/scratch/molevo/dts/os/dna/hic/omnic/P21656_102/02-FASTQ/211020_A00621_0515_BHLC72DSX2/P21656_102_S7_L002_R1_001.fastq.gz
R2=/scratch/molevo/dts/os/dna/hic/omnic/P21656_102/02-FASTQ/211020_A00621_0515_BHLC72DSX2/P21656_102_S7_L002_R2_001.fastq.gz
MYBIN=/home/user/destanovic/bin/


${MYBIN}/chromap/chromap -i -r ${ASSEM} -o ${PREFIX}.index

awk -f ${MYBIN}/generate-assembly-file-from-fasta.awk ${ASSEM} > ${PREFIX}.assembly

module unload python
module load python2
python ${MYBIN}/assembly-from-fasta.py -c ${ASSEM} py_${PREFIX}
module unload python
module python3

for qual in {2,10}
do
${MYBIN}/chromap/chromap -t ${THREAD} --preset hic -x ${PREFIX}.index -r ${ASSEM} -1 ${R1} -2 ${R2} -q ${qual} -o ${PREFIX}_${qual}.pairs
grep -v '#' ${PREFIX}_${qual}.pairs | \
        awk '{{if($6!="+") $6=16; else $6=0; if($7!="+") $7=16; else $7=0}} \
            $2<=$4{{print $6, $2, $3, 0, $7, $4, $5, 1, "1 - - 1  - - -" }} \
            $4<$2{{print $7, $4, $5, 0, $6, $2, $3, 1, "1 - - 1  - - -" }}' > ${PREFIX}_${qual}.longp
bash ${MYBIN}/3d-dna/visualize/run-assembly-visualizer.sh -p false ${PREFIX}.assembly ${PREFIX}_${qual}.longp || true
mv ${PREFIX}.hic ${PREFIX}_hic_${qual}.hic
done
