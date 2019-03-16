#!/bin/bash
set -e
set -u
set -o pipefail

#$1 is filepath to pipeline output folder

outdir=${1}
pidcutoff=${2}
qcovcutoff=${3}

> ${1}/blast/allalignments.tsv
echo -e "qseqid\tsseqid\tpident\tevalue\tbitscore\tqcovhsp\tqlen\tslen" >> ${1}/blast/allalignments.tsv 
#reformat blast output file                                                                                                                
samples=($(cut -f1 ${outdir}/blastdbfilepaths.tsv | sort -V | uniq)) #$subject sequence blast database file; format: sample \t filepath to database

for sample in ${samples[@]}; do
    mkdir -p ${outdir}/blast/${sample}
    cat ${outdir}/blast/${sample}_alignments.tsv | awk -v "var1=${pidcutoff}" -v "var2=${qcovcutoff}" '$3 >= var1 && $6 >= var2' | tee ${outdir}/blast/${sample}/alignments.tsv >> ${outdir}/blast/allalignments.tsv #percent identity >=40% ; hsp query coverage >= 80%
    rm ${outdir}/blast/${sample}_alignments.tsv
done



#OLD CODE

#don't want to output included at this stage since they may not be included when considering reciprocal best hits only
#> ${1}/included.txt
#| tee ${1}/blast/${sample}/alignments.tsv | if [ $(wc -l) -gt 0 ]; then echo "${sample}" >> ${1}/included.txt; fi



#| awk '$3 > 40 && $8 >80' | awk '{OFS="\t"; if($6 < $7) {print $0"\t""+"} else {print $0"\t""-"}}'
# awk '{OFS="\t"; if($9 < $10) {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; print $0"\t""+"} else {qcovprop=($13/100); qcovhspprop=($14/100); $13=qcovprop; $14=qcovhspprop; sstart=$10; send=$9; $9=sstart; $10=send; print $0"\t""-"}}' 

     
