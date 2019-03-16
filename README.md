metamorth is a command-line tool for identifying [orthologs](https://en.wikipedia.org/wiki/Sequence_homology#Orthology) amongst pairs of bacterial genomes, and calculating re-arrangement distances, according to the pairwise ordering of orthologs. The [reciprocal best hit approach](https://www.ncbi.nlm.nih.gov/pubmed/23160176) is used to identify hypothetical orthologs between genomes, whilst breakpoint distance is calculated as a [re-arrangement distance metric](https://cse.sc.edu/~jtang/mage.pdf). metamorth performs an all-vs-all comparison (i.e. orthologs/breakpoint distances are determined for each pair of genomes).


# Table of contents

* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Input](#Input)
* [Quick start](#Quick-start)
* [Background and methods](#Background-and-methods)
* [Options and usage](#Options-and-usage)
* [Output files](#Output-files)
* [License](#License)



# Introduction

The reciprocal best hit approach is a widely used method to identify candidate orthologs among pairs of genomes. Briefly, BLAST alignments are conducted between the sets of proteins encoded by a pair of genomes; proteins which produce mutual best hit alignments in both blast directions (genome A vs genome B and genome B vs genome A) are considered to be hypothetical orthologs. Since all-vs-all BLAST alignments can be time-consuming, the reciprocal best hit component of metamorth is intended for small-to-medium-sized datasets (e.g. plasmids; small collections of bacterial genomes). If you wish to analyse large datasets, [faster algorithms](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101850) can be used to infer orthologs; in this case, metamorth can still be used to determine breakpoint distances. 

As input, metamorth can be provided with one of the following:

* A Genbank file ([full format](https://widdowquinn.github.io/2018-03-06-ibioic/01-introduction/02-annotation.html), i.e. with gene locus information) comprising nucleotide sequences of the genomes that are to be compared. For details, see the [Input](#Input) section.
* A file containing orthologs (determined using the user's preferred method), formatted as described in [Input](#Input).


# Requirements

* Linux or MacOS (with the [Bash shell](https://en.wikibooks.org/wiki/Bash_Shell_Scripting#What_is_Bash?), which is the default shell on MacOS and many Linux distributions)
* [Python](https://www.python.org/) 2.7 or Python 3
* [SeqKit](https://github.com/shenwei356/seqkit)
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (`blastp`)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [R](https://www.r-project.org/) 3.3.1 or later with the following packages installed:
    * [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html); [gsubfn](https://cran.r-project.org/web/packages/gsubfn/index.html); [foreach](https://cran.r-project.org/web/packages/foreach/index.html); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html); [data.table](https://cran.r-project.org/web/packages/data.table/index.html); [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)<br>

Run the following code in R to install the required R packages:<br>
```bash
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

install.packages("devtools",repo='https://cloud.r-project.org/')
library(devtools)
devtools::install_github("ggrothendieck/gsubfn")

install.packages("tidyr",repo='https://cloud.r-project.org/')
install.packages("foreach",repo='https://cloud.r-project.org/')
install.packages("doParallel",repo='https://cloud.r-project.org/')
install.packages("data.table",repo='https://cloud.r-project.org/')
```

# Installation

```bash
git clone https://github.com/AlexOrlek/metamorth.git
cd metamorth
```
You should find the metamorth.py executable script within the repository directory. If you add the path of this directory to your [$PATH variable](https://www.computerhope.com/issues/ch001647.htm), then metamorth can be run by calling `metamorth.py [`*`arguments...`*`]` from any directory location. Note also that metamorth expects the tools listed in [Requirements](#Requirements) to be available in your $PATH variable.

# Input

### Genbank file input

If a genome has already been uploaded to Genbank, a full Genbank file can be downloaded with the following code, after having installed [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/), replacing 'queryname' with an NCBI accession id:<br>
`esearch -db nuccore -query queryname | efetch -format gb -style withparts > sequence.gb`<br>
Alternatively, for newly sequenced genomes, annotation software can be run, with the output format specified as Genbank. For example, [prokka](https://github.com/tseemann/prokka) annotation software provides Genbank output format as an option.


### Best hit alignments / ortholog input

If orthologs have already been determined, these can be provided directly to metamorth, and breakpoint distances can then be calculated. The ortholog input file must be formatted as in the example file contained in the /example folder. There should be 2 columns of data representing reciprocal hits in the format: genome_id|protein_id|strand|nucleotide_position. If orthologs have not been determined, but best hits have been determined, then these hits can be provided, following the same formatting rules; metamorth ensures only reciprocal best hits (candidate orthologs) are retained. A header row at the top of the input file must be present, providing a name for each column (there is no rule about how the columns should be named).


# Quick start

The complete pipeline (protein extraction, ortholog determination, breakpoint distance calculation) can be run as follows, by providing a Genbank file using the `-s` flag:

`metamorth.py -s sequences.gb -o output-directory --breakpoint`



# Options and usage

`metamorth.py --help` produces a summary of all the options.


# Background and methods

Further general information about ortholog detection can be found in a recent [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5674930/) by Nichio _et al_. A brief outline of the steps of the complete metamorth pipeline is given below:

1. Translated protein sequences are extracted from a Genbank file containing annotated nucleotide sequences of genomes. In addition, the nucleotide position of the corresponding coding sequence, and the strand (+ve/-ve) of the translated protein sequence are extracted.
2. BLAST is conducted between the sets of protein sequences from each pair of genomes, in both BLAST directions (i.e. genome A proteins vs genome B proteins and genome B proteins vs genome A proteins). For each protein, only the best hit across all proteins in the other genome is retained, and only if this hit satisfies identity and coverage thresholds (default: at least 40% protein sequence identity and at least 80% query coverage).
3. The remaining BLAST hits are further filtered, retaining hits only if they are reciprocal best hits.
4. The reciprocal best hits (candidate orthologs) are ordered by nucleotide start position, generating a [signed permutation](http://rosalind.info/glossary/signed-permutation/) for each genome in a given pairwise comparison (i.e. a sequence of numbers, associated with +/- signs, representing the ordered orthologs and their strand). From these signed permutations, a [breakpoint distance](https://www.liebertpub.com/doi/abs/10.1089/cmb.1998.5.555) is calculated and expressed as number of breakpoints / number of shared orthologs.



# Output files

The below table shows the outputs from running the complete pipeline (protein extraction, ortholog determination, breakpoint distance calculation).

File/Directory               | Description
---------------------------- | -------------------------------------------------------------------------------------------------
fastafilepaths.tsv	     | genome names and corresponding FASTA file paths
blastdbfilepaths.tsv	     | genome names and corresponding BLAST database file paths
fastas/                      | directory containing protein FASTA files (and corresponding BLAST databases), derived from the input Genbank file (one FASTA file per genome)
fastas/allfastas.fasta       | concatenated protein FASTA files
fastas/sequencestats.tsv     | file containing statistics on number of proteins extracted from coding sequences of each genome
blast/                       | directory containing tsv files of best hit BLAST alignments for each genome
blast/allalignments.tsv      | best BLAST hits for all pairwise combinations
blast/allalignments_RBH.tsv  | reciprocal best BLAST hits (candidate orthologs) for all pairwise combinations
output/                      | directory containing the breakpointdistance.tsv file (pairwise breakpoint distances)
included.txt                 | file showing names of genomes included in the analysis (any genome sharing at least 1 reciprocal best hit with another genome)



# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
