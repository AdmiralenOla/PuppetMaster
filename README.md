# PuppetMaster
Tools for manipulating sequencing data, multiple sequence alignments and phylogenetic trees

# VARIABLE SITES EXTRACTION

The aim of many microbial typing pipelines is to compare isolates at the SNP level. After multiple sequence alignment,
phylogenetic inferences can be made by comparing SNPs across isolates. However, since SNPs are usually called with respect
to a reference genome, one typically ends up with many redundant SNP sites when comparing isolates towards each other.

In parsimony methods (as opposed to likelihood or distance methods), constant sites (i.e. sites where the isolates being
compared do not differ) can safely be excluded, as they are not informative for tree inference.

Similarly, sites where just one isolate is polymorphous, although variable, may not be of interest since they do not
contribute to tree discrimination. Finally, sites where all four bases are represented are also non-informative and
should be trimmed.

Furthermore, multiple sequence alignments are typically cluttered with gap characters "-" and ambigous characters "N".
These sites are not valuable for phylogenetic inference, and should be trimmed.

This script allows the user to trim away N-containing columns, gapped columns and non-variable/non-informative sites
from a multiple sequence alignment FASTA file.

Author: Ola Brynildsrud
