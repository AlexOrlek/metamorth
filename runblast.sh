#!/bin/bash
set -e
set -u
set -o pipefail

mkdir -p ${1}/blast
mkdir -p ${1}/output

query=${1}/fastas/allfastas.fasta #query fasta file (all samples)
blastdbfilepaths=${1}/blastdbfilepaths.tsv #blastdbfilepaths=${1}/blastdbfilepaths_TEST.tsv
threads=${3}
evalue=${2} #default: '1e-6' #see Moreno-Hagelsieb 2008 Choosing BLAST options...
sort -k1,1V -o ${blastdbfilepaths} ${blastdbfilepaths}


while IFS=$'\t' read -r -a line
do
    sample="${line[0]}"
    database="${line[1]}"
    blastoutput="${1}/blast/${sample}_alignments.tsv"
    cat ${query} | seqkit grep -r -p ^${sample} -v  | blastp -db ${database} -out ${blastoutput} -outfmt '6 qseqid sseqid pident qcovhsp qlen slen' -task 'blastp' -num_threads ${threads} -max_target_seqs '1' -max_hsps '1' -evalue ${evalue}
done < ${blastdbfilepaths}


#OLD CODE
#not sure how to redirect stdout to pipe - just reformat at later stage
# | awk '$3 > 40 && $8 >80' | awk '{OFS="\t"; if($6 < $7) {print $0"\t""+"} else {print $0"\t""-"}}'


#old outfmt - don't need so much output
#-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen' -task 'blastp' -num_threads ${threads} -max_target_seqs '1' -max_hsps '1' -evalue ${evalue} | awk '$3 > 40 && $14 >80' 
#-word_size ${wordsize} -culling_limit '5'
#-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen'
#-evalue ${evalue}


#| filterbydifference.py
