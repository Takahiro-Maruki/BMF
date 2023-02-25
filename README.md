# BMF
February 25, 2023

Author: Takahiro Maruki

Software package of the Bayesian mutation finder (BMF). 

**This pacakge can be downloaded by clicking "Code" at this web site.**

This software is for identifying single nucleotide mutations from high-throughput sequencing data from mutation accumulation (MA) experiments.  It calls genotypes in MA lines using BGC (Maruki and Lynch 2017) and filters individual genotype calls using a binomial test examining the goodness of fit between observed nucleotide read counts and their binomial expectation.

**Input file**.  The input is a tab-delimited text file and can be made from a pro file of nucleotide read counts using GFE (Maruki and Lynch 2015; GFE version 3.0 included in the Programs directory in this software package) in the `c’ mode (available from GFE version 2.0).  The meanings of the first twelve columns are: 1) scaffold (chromosome) identifier; 2) site identifier (coordinate); 3) nucleotide of the reference sequence; 4, 5) nucleotides of the major and minor alleles, respectively (1: A, 2: C, 3: G, 4: T); 6) depth of coverage in the population sample (sum of the coverage over the individuals; 7) error rate estimate; 8, 9, 10) ML estimates of the frequencies of major homozygotes, heterozygotes, and minor homozygotes, respectively; 11) likelihood-ratio test statistic for polymorphism; 12) likelihood-ratio test statistic for deviation from Hardy-Weinberg equilibrium.  Thereafter, the counts of nucleotide reads A, C, G, T separated by slashes are shown for each individual in each of the columns.  The likelihood-ratio test statistics for polymorphism and HWE-deviation are expected to asymptotically correspond to chi-squared distributions with two and one degrees of freedom, respectively.

**Output file**.  The output file is also a tab-delimited text file.  The meanings of the first eight columns are: 1) scaffold (chromosome) identifier; 2) site identifier (coordinate); 3) nucleotide of the reference sequence; 4) error rate estimate; 5) depth of coverage in the population sample (sum of the coverage over the individuals); 6) inferred ancestral genotype; 7) number of genotype calls; 8) minor-allele count.  Thereafter, the genotype call is shown for each individual in each of the columns.

**Prerequisite**

The installation of the following is necessary.  gcc and Perl are generally pre-installed in a Linux server, where analysis of high-throughput sequencing data is generally carried out.

1. gcc (https://gcc.gnu.org/)

2. Perl (https://www.perl.org/)

3. Samtools (Li et al. 2009) (http://www.htslib.org/)

**Installation**

Compile the Perl script and C++ programs in the Programs directory.

1. Ext_Ref_Nuc.pl <br />
chmod u+x Ext_Ref_Nuc.pl

2. mpileup2pro.cpp <br />
g++ -o mpileup2pro mpileup2pro.cpp -lm

3. HGC.cpp <br />
g++ -o HGC HGC.cpp -lm

4. Rem_Multi_Allelic.cpp <br />
g++ -o Rem_Multi_Allelic Rem_Multi_Allelic.cpp -lm

5. GFE_v3.0.cpp <br />
g++ -o GFE_v3.0 GFE_v3.0.cpp -lm

6. bmf.cpp <br />
g++ -o bmf bmf.cpp -lm

**Usage instructions**

In the following, the script and programs are assumed to be in $PATH.  If thay are not in $PATH, please specify the path (location) of each script/program to run it.

1. Make an mpileup file from BAM files using the mpileup function of Samtools (Li et al. 2009). <br />
samtools mpileup -b BamList.txt -f Reference.fasta -o Out.mpileup

2. Extract the nucleotide at every position of the reference sequence using Ext_Ref_Nuc.pl. <br />
perl Ext_Ref_Nuc.pl Reference.fasta RefNuc.txt

3. Make a pro file of nucleotide read counts from the mpileup file using mpileup2pro. <br />
mpileup2pro -ref RefNuc.txt -id IDs.txt -mp Out.mpileup -out Out.pro
- The order of individuals listed in IDs.txt needs to be consistent with that in the mpileup file.

4. Run HGC (Maruki and Lynch 2017) to identify tri- and tetra-allelic sites, setting the minimum required coverage to call individual genotypes at six. <br />
HGC -in Out.pro -min_cov 6 -out Out_mc6_HGC.txt

5. Set the depth of coverage of all individuals to zero in the pro file at tri- and tetra-allelic sites using Rem_Multi_Allelic. <br />
awk -v OFS='\t' '{if ($1 == "scaffold" || $4 != "NA" && $4 > 2) print $1, $2}' Out_mc6_HGC.txt > List_MA_Out_mc6_HGC.txt <br />
Rem_Multi_Allelic -pf Out.pro -mf List_MA_Out_mc6_HGC.txt -out MAR_Out.pro

6. Run GFE_v3.0 (Maruki and Lynch 2015) in the c mode. <br />
GFE_v3.0 -in MAR_Out.pro -mode c -out Out_c_GFEv3.0.txt

7. Run bmf, setting the minimum required coverage to call individual genotypes at eight and outputting reference nucleotides. <br />
bmf -in Out_c_GFEv3.0.txt -rn 1 -min_cov 8 -out Out_rn1_mc8_bmf.txt
- The `-in`, and `-out` options specify the input and output file names, respectively. <br />
- The reference nucleotides can be shown in the output by setting the `-rn' option at one (1).
- The minimum and maximum coverage for calling individual genotypes can be specified by the `-min_cov` and `-max_cov`, respectively.  The default values for the minimum and maximum coverage are 1 and 2,000,000,000, respectively. <br />
- The critical values for the heterozygous and homozygous cumulative binomial probabilities in the binomial test can be specified by `-cv_het` and `-cv_hom` options, respectively.  The default values for the heterozygous and homozygous cumulative binomial probabilities are 0.025 and 0.05, respectively.
- Results at all sites in the reference sequence can be shown in the output file by setting the `-as` option at one (1).  The input file also needs to contain all sites in the reference sequence in this case, which can be done by running GFE_v3.0 in the `c` mode and setting the `-as` option at one (1). <br />
- A usage help message explaining these options can be shown by typing the following command: <br />
bmf -h

**Reference**

If you use this software, please cite the following paper:

Maruki, T, Ozere, A, and Cristescu, M. E., Genome-wide identification of single nucleotide mutations from time-series mutation accumulation data. In prep.

**Copyright Notice**

This program is freely available; and can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
For a copy of the GNU General Public License write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

**Contact**

If you have difficulty using this software, please send the following information to Takahiro Maruki (takahiro.maruki@mcgill.ca):

1. Brief explanation of the problem.

2. Command used.

3. Part of the input file.

4. Part of the output file.
