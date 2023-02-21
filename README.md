# BMF
February 19, 2023

Software package of the Bayesian mutation finder (BMF). 

**This pacakge can be downloaded by clicking "Code" at this web site.**

This software is for identifying mutations from high-throughput sequencing data from mutation accumulation (MA) experiments.  It calls genotypes in MA lines using BGC (Maruki and Lynch 2017) and filters individual genotype calls using a binomial test examining the goodness of fit between observed nucleotide read counts and their binomial expectation.

**Input file**.  The input is a tab-delimited text file and can be made from a pro file of nucleotide read counts using GFE (Maruki and Lynch 2015; GFE version 3.0 included in the Programs directory in this software package) in the `c’ mode (available from GFE version 2.0).  The meanings of the first twelve columns are: 1) scaffold (chromosome) identifier; 2) site identifier (coordinate); 3) nucleotide of the reference sequence; 4, 5) nucleotides of the major and minor alleles, respectively (1: A, 2: C, 3: G, 4: T); 6) depth of coverage in the population sample (sum of the coverage over the individuals; 7) error rate estimate; 8, 9, 10) ML estimates of the frequencies of major homozygotes, heterozygotes, and minor homozygotes, respectively; 11) likelihood-ratio test statistic for polymorphism; 12) likelihood-ratio test statistic for deviation from Hardy-Weinberg equilibrium.  Thereafter, the counts of nucleotide reads A, C, G, T separated by slashes are shown for each individual in each of the columns.  The likelihood-ratio test statistics for polymorphism and HWE-deviation are expected to asymptotically correspond to chi-squared distributions with two and one degrees of freedom, respectively.

**Output file**.  The output file is also a tab-delimited text file.  The meanings of the first eight columns are: 1) scaffold (chromosome) identifier; 2) site identifier (coordinate); 3) nucleotide of the reference sequence; 4) error rate estimate; 5) depth of coverage in the population sample (sum of the coverage over the individuals); 6) inferred ancestral genotype; 7) number of genotype calls; 8) minor-allele count.  Thereafter, the genotype call is shown for each individual in each of the columns.

**Reference**

If you use this software, please cite the following paper:

Maruki, T, Ozere, A, and Cristescu, M. E., Genome-wide identification of single nucleotide mutations from time-series mutation accumulation data. In prep.

