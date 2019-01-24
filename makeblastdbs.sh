#!/bin/bash
set -e
set -u
set -o pipefail

#arg[1] is outputpath; arg[2] is threads; arg[3] is sourcedir

#$1 is fastafilepaths.tsv

#concatenate fasta files
> ${1}/fastas/allfastas.fasta
cat ${1}/fastafilepaths.tsv | cut -f2 | xargs cat >> ${1}/fastas/allfastas.fasta
    
#make blast database of each individual fasta file
cat ${1}/fastafilepaths.tsv | cut -f2 | python ${3}/removeextension.py | parallel -k -j ${2} "makeblastdb -dbtype prot -in {}.fasta -out {}_db"
