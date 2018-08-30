PopPoly, version 2.2, Ehsan Motazedi (ehsan.motazedi@gmail.com), last modified: Aug 15 2018

Introduction
============

**PopPoly** \(executable *main.py*\) is a tool for haplotype estimation in polyploid F1 populations. PopPoly uses all of the population reads to estimate the maximum likelihood parental haplotypes, from which offspring haplotypes are selected using the MEC criterion.

The algorithm is based on a Bayesian framework to choose the most likely parental haplotypes using sequence reads, and applies Mendelian inheritance rules, assuming no recombination between the SNPs, to weight each candidate phasing for the parents and to later choose the offspring phasings from the parental haplotypes. It can also re-assign the dosages for all or missing SNPs using the Mendelian rules and sequence reads. The name of the parents in the F1 population must be provided as input parameter, with the possibility to specify just one parent in self-fertilizing populations. Otherwise the first two samples registered in the BAM header will be counted as parents. Also, the VCF file must contain only SNPs, that is to say no complex variants and indels. This can be achieved by running the *break_vcf.sh* script on the original VCF file.

PopPoly is written in *Python* and is compatible with both Python 2 and Python 3. There is also a perl wrapper available *PopPoly_wrapped.pl* to directly use input files that contain multiple contigs, in which case one has to specify the name of the desired contig for haplotyping through the wrapper option. The wrapper also allows to filter SNPs based upon their genotype missing rate in the population and can be run using the original VCF file as it calls break_vcf.sh internally (see PopPoly_wrapped.pl -h).

It is also possible to exclude some of the individuals from analysis, for example if several replicates of the same sample are present under different names and one wishes to use only one of them. The space separated names of the to be excluded individuals must be stord in the shell environment in a variable called *exclude* before running PopPoly, as follows:

export exclude="sample_1_replciate_2 sample_2_replicate_1 sample_3_replicate_2 sample_3_replicate_3" 

Note that PopPoly relies on *pysam* and *numpy*, which must have been properly installed before PopPoly is called. Some users have reported problems with *numpy undefined symbol*, a discussion of which can be found [here] (https://stackoverflow.com/questions/36190757/numpy-undefined-symbol-pyfpe-jbuf).

Input parameters:
=====================

For the input parameters, see the help of the software:

./main.py --help, -h

***Example:***


./main.py bamfile vcffile output_dirname --kappa 0.85 -a 0.8 --rho 0.2 --impute

Citation:
=====================
To cite PopPoly, please refer to *Family-based haplotype estimation and allele dosage correction for polyploids using short sequence reads*, Ehsan Motazedi, Chris Maliepaard, Richard Finkers and Dick de Ridder, 2018, bioRxiv 318196; *doi: https://doi.org/10.1101/318196*

Download:
=====================
The software is available at gitlab, Wageningen UR:

git clone git@git.wageningenur.nl:motaz001/PopPoly.git —recursive

Copyright and License:
=====================
Copyright (c) 2017, Ehsan Motazedi, Wageningen UR.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files \(the "Software"\), to deal in the Software without restriction, including without limitation the rights to use, copy, share and modify the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. This permission is only valid for academic and personal use, i.e. not for commercial purposes.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
