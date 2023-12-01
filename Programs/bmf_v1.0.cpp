/* Updated on 11/11/23 */

/*

Program bmf_v1.0.cpp to estimate the genotype of each
MA line form pro files of individual high-throughput
sequencing data of mutation-accumulation experiments,
incorporating pre-determined genotype-frequency estimates
by the ML genotype-frequency estimator.  Genotype probabilities
for each individual can be shown.  The output at all sites corresponding
to the input at all sites can be printed out by setting the `-as' option
at one.  Reference nucleotides can be printed out by setting the `-rn'
option at one.  The ancestral genotype is inferred from the genotype
frequency estimates.  Genotype calls are refined using a binomial test
assuming the inferred ancestral genotype in the null hypothesis.  In this
version 1.0, the error rate estimate is included in the equation for the
binomial test with a heterozygous genotype.

*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/types.h>
#include <vector>
using namespace std;

// Function declaration
double binomial(int n, int x, double p);

// Function returning a binomial probability
double binomial(int n, int x, double p) {
	// Local variable declarations
	double ln_p;	// natural logarithm of the binomial probability
	double prob;	// binomial probability
	int ag, bg, cg;	// counters

	ln_p = 0.0;
	for (ag=n; ag>=1; ag--) {
		ln_p = ln_p + log(ag);
	}
	for (bg=x; bg>=1; bg--) {
		ln_p = ln_p - log(bg);
	}
	for ( cg=(n-x); cg>=1; cg-- ) {
		ln_p = ln_p - log(cg);
	}
	ln_p = ln_p + x*log(p) + (n - x)*log( (double)1.0 - p );
	prob = exp(ln_p);
	return prob;
}

int main(int argc, char *argv[]) {

	// Default values of the options
	const char* in_file_name = {"In_bmf.txt"};
	const char* out_file_name = {"Out_bmf.txt"};
	int min_cov = 1;	// minimum coverage cut-off required to call a genotype
	int max_cov = 2000000000;	// maximum coverage cut-off allowed to call a genotype
	int gp = 0;	// indicator to specify whether genotype probabilities are shown (0: do not show, 1: show)
	int as = 0;	// indicator to specify whether the output at all sites is shown (0: do not show, 1: show)
	int rn = 1;	// indicator to specify whether the output includes reference nucleotides (0: do not include, 1: include)
	double cv_het = 0.025;	// critical value for heterozygous cumulative binomial probability
	double cv_hom = 0.05;  // critical value for the homozygous cumulative binomial probability
	int print_help = 0;

	int argz = 1;	// argument counter

	// Read the specified setting
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-min_cov") == 0) {
			sscanf(argv[++argz], "%d", &min_cov);
		} else if (strcmp(argv[argz], "-max_cov") == 0) {
                        sscanf(argv[++argz], "%d", &max_cov);
		} else if (strcmp(argv[argz], "-gp") == 0) {
			sscanf(argv[++argz], "%d", &gp);
		} else if (strcmp(argv[argz], "-as") == 0) {
                        sscanf(argv[++argz], "%d", &as);
		} else if (strcmp(argv[argz], "-rn") == 0) {
                        sscanf(argv[++argz], "%d", &rn);
		} else if (strcmp(argv[argz], "-cv_het") == 0) {
			sscanf(argv[++argz], "%lf", &cv_het);
		} else if (strcmp(argv[argz], "-cv_hom") == 0) {
			sscanf(argv[++argz], "%lf", &cv_hom);
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) {	// print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options: -h: print the usage message\n");
		fprintf(stderr, "	-in <s>: specify the input file name\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-min_cov <d>: specify the minimum coverage cut-off required to call a genotype\n");
		fprintf(stderr, "       -max_cov <d>: specify the maximum coverage cut-off allowed to call a genotype\n");
		fprintf(stderr, "	-gp <d>: specify whether genotype probabilities are shown\n");
		fprintf(stderr, "       -as <d>: specify whether the output at all sites is shown\n");
		fprintf(stderr, "       -rn <d>: specify whether the output includes reference nucleotides\n");
		fprintf(stderr, "       -cv_het <f>: specify the critical value for the heterozygous cumulative binomial probability\n");
		fprintf(stderr, "       -cv_hom <f>: specify the critical value for the homozygous cumulative binomial probability\n");
		exit(1);
	}

	if (min_cov > 1) {
		printf("min_cov: %d\n", min_cov);
	}
	if (max_cov < 2000000000) {
		printf("max_cov: %d\n", max_cov);
	}

	string line;	// String buffer

	ifstream inputFile(in_file_name);	// Try to open the input file
	if ( !inputFile.is_open() ) {	// Exit on failure
		printf("Cannot open %s for reading.\n", in_file_name);
		exit(1);
	}

	// Read the header
	string h_scaf, h_site, h_ref_nuc, h_n1, h_n2, h_pop_cov, h_best_err, h_best_P, h_best_H, h_best_Q, h_pol_llstat, h_HWE_llstat;	// Stores header names
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_scaf >> h_site >> h_ref_nuc >> h_n1 >> h_n2 >> h_pop_cov >> h_best_err >> h_best_P >> h_best_H >> h_best_Q >> h_pol_llstat >> h_HWE_llstat;
	string str;	// Temporarily stores each individual ID
	vector <string> id_ind;      // Stores individual IDs.
	id_ind.clear();
	while (true) {
		ss >> str;
		id_ind.push_back(str);
		if ( ss.eof() ) {
			break;
		}
	}
	int nsample = id_ind.size();
	printf("%d individuals to be analyzed\n", nsample);

  	int ig, gg, jg, kg, lg, mg, size_quartet, digit, count_nuc;
	string scaffold, ref_nuc;
	int site;
	double best_error, best_error13, best_geno_freq[4], best_p;
  	string s_n1, s_n2, s_best_error, s_best_geno_freq[4], s_pol_llstat, s_HWE_llstat, quartet[nsample+1];
  	int ind_coverage[nsample+1], num_digit[5], read_count[nsample+1][5];
	int pop_coverage, n1, n2;
	string major_allele, minor_allele, best_geno[nsample+1];
	int ml_ind_read[nsample+1][4], count_geno[4], tot_count_geno, major_allele_count, minor_allele_count, id_best_geno, ml_geno;
	double sum_weights, maxll, sub_maxll, prob_nuc[4][4], prob_data[nsample+1][4], llg[4];
	string anc;	// ancestral genotype
	double prob_ind_het, prob_ind_hom;      // cumulative probabilities of heterozygous and homozygous individual genotypes
	int dg;	// depth counter
	int f_count_geno[4];	// genotype call counts excluding those containing mutations
	int f_major_allele_count, f_minor_allele_count;	// filtered allele counts
	int num_mut_lines;	// number of MA lines with mutations
	double base_prob_ind_hom;	// baseline probability of a homozygous individual genotype
	string f_anc;	// ancestral genotype inferred from genotype calls in MA lines without mutation

	// Point to the output file
	FILE *outstream;

  	// Open the output file
  	outstream = fopen(out_file_name, "w");

  	// Print the field names
	if (rn == 0) {
		fprintf(outstream, "scaffold\tposition\tbest_error\tpop_coverage\tanc\tn_gc\tmac\t");
		// printf("scaffold\tposition\tbest_error\tpop_coverage\tanc\tn_gc\tmac\t");
	} else if (rn == 1) {
		fprintf(outstream, "scaffold\tposition\tref_nuc\tbest_error\tpop_coverage\tanc\tn_gc\tmac\t");
		// printf("scaffold\tposition\tref_nuc\tbest_error\tpop_coverage\tanc\tn_gc\tmac\t");
	}
	if (gp == 0) {
		for(ig=0;ig<nsample-1;ig++) {
			fprintf(outstream, "%s\t", id_ind[ig].c_str());
			// printf("%s\t", id_ind[ig].c_str());
		}
		fprintf(outstream, "%s\n", id_ind[ig].c_str());
		// printf("%s\n", id_ind[ig].c_str());
	} else if (gp == 1) {
		for(ig=0;ig<nsample-1;ig++) {
			fprintf(outstream, "%s\t%s GP\t", id_ind[ig].c_str(), id_ind[ig].c_str());
        		// printf("%s\t%s GP\t", id_ind[ig].c_str(), id_ind[ig].c_str());
		}
		fprintf(outstream, "%s\t%s GP\n", id_ind[ig].c_str(), id_ind[ig].c_str());
		// printf("%s\t%s GP\n", id_ind[ig].c_str(), id_ind[ig].c_str());
	}
  	// Read the main data
  	while ( getline(inputFile, line) ) {
  		istringstream ss(line);
		ss >> scaffold >> site >> ref_nuc >> s_n1 >> s_n2 >> pop_coverage >> s_best_error >> s_best_geno_freq[1] >> s_best_geno_freq[2] >> s_best_geno_freq[3] >> s_pol_llstat >> s_HWE_llstat;
		if (s_best_geno_freq[1] != "NA") {
			n1 = atoi(s_n1.c_str());
			n2 = atoi(s_n2.c_str());
			best_error = atof(s_best_error.c_str());
			best_geno_freq[1] = atof(s_best_geno_freq[1].c_str());
			best_geno_freq[2] = atof(s_best_geno_freq[2].c_str());
			best_geno_freq[3] = atof(s_best_geno_freq[3].c_str());
			best_p = best_geno_freq[1] + best_geno_freq[2]/2.0;
			for (ig = 1; ig <= nsample; ig++) {
				ss >> quartet[ig];
			}
			for (ig = 1; ig <= nsample; ig++) {
				size_quartet = quartet[ig].size();
				// printf("quartet size: %d\n", size_quartet);
				jg = 1;
				digit = 0;
				ind_coverage[ig] = 0;
				for (kg = 0; kg < size_quartet; kg++) {
					if ( quartet[ig].at(kg) == '/') {
						// cout << quartet[ig].at(kg) << endl;
						mg = pow(10,digit-1);
						// printf("digit: %d\tmg: %d\n", digit, mg);
						// cout << digit << endl;
						count_nuc = 0;
						lg = 0;
						while (lg <= digit-1) {
							// cout << "Entered the for loop" << endl;
							// cout <<  << endl;
							count_nuc = count_nuc + mg*num_digit[lg];
							mg = mg/10;
							num_digit[lg] = 0;
							lg = lg + 1;
						}
						// printf("jg: %d\n", jg);
						read_count[ig][jg] = count_nuc;
						ind_coverage[ig] = ind_coverage[ig] + count_nuc;
						jg = jg + 1;
						count_nuc = 0;
						digit = 0;
					} else {
						// cout << quartet[ig].at(kg) << endl;
						num_digit[digit] = quartet[ig].at(kg) - '0';
						// cout << num_digit[digit] << endl;
						digit = digit + 1;
					}
					if (kg == size_quartet - 1) {
						mg = pow(10,digit-1);
						// printf("digit: %d\tmg: %d\n", digit, mg);
						count_nuc = 0;
						lg = 0;
						while (lg <= digit-1) {
							// cout << lg << endl;
							count_nuc = count_nuc + mg*num_digit[lg];
							mg = mg/10;
							num_digit[lg] = 0;
							lg = lg + 1;
						}
						// printf("jg: %d\n", jg);
						read_count[ig][jg] = count_nuc;
						ind_coverage[ig] = ind_coverage[ig] + count_nuc;
						jg = 1;
					}
				}
			}
			// cout << scaffold << site << ref_nuc.c_str() << quartet[1].c_str() << endl;
			// read_count[1][1] = quartet[1].at(0) - '0';
			// cout << read_count[1][4] << endl;
			/*
			for (ig = 1; ig <= nsample; ig++) {
				printf("%s\t%d\t%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str(), ig, id_ind[ig].c_str());
				fprintf(outstream, "%s\t%d\t%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str(), ig, id_ind[ig].c_str());
				for (jg = 1; jg <= 4; jg++) {
					printf("%d\t", read_count[ig][jg]);
					fprintf(outstream, "%d\t", read_count[ig][jg]);
				}
				printf("%d\n", ind_coverage[ig]);
				fprintf(outstream, "%d\n", ind_coverage[ig]);
			}
			*/
			if (best_p < 1.0) {
				// ML estimation starts here
        			// Assign nucleotides to the candidate alleles at the site.
				if (n1 == 1) {
					major_allele = "A";
				} else if (n1 == 2) {
					major_allele = "C";
				} else if (n1 == 3) {
					major_allele = "G";
				} else if (n1 == 4) {
					major_allele = "T";
				}
				if (n2 == 1) {
					minor_allele = "A";
				} else if (n2 == 2) {
					minor_allele = "C";
				} else if (n2 == 3) {
					minor_allele = "G";
				} else if (n2 == 4) {
					minor_allele = "T";
				}
				// Infer the ancestral genotype using the genotype frequency estimates
				if (best_geno_freq[1] > best_geno_freq[2]) {
					anc = major_allele + major_allele;
				} else {
					anc = major_allele + minor_allele;
				}
				// The default individual genotypes are NA
				for (ig=1;ig<=nsample;ig++) {
					best_geno[ig] = "NA";
					prob_data[ig][1] = 0.0;
					prob_data[ig][2] = 0.0;
					prob_data[ig][3] = 0.0;
					// printf("%s\n", best_geno[ig].c_str());
				}
				// Carry out the ML estimation
				// printf("Just before the ML estimation\n");
				count_geno[1] = 0;
  				count_geno[2] = 0;
  				count_geno[3] = 0;
				major_allele_count = 0;
				minor_allele_count = 0;
  				tot_count_geno = 0;
				f_count_geno[1] = 0;
                                f_count_geno[2] = 0;
                                f_count_geno[3] = 0;
                                f_major_allele_count = 0;
                                f_minor_allele_count = 0;
                                num_mut_lines = 0;
				// Calculate the probability of a nucleotide read given the genotype of the individual
				best_error13 = best_error/3.0;
				prob_nuc[1][1] = 1.0 - best_error;
				prob_nuc[1][2] = best_error13;
				prob_nuc[1][3] = best_error13;
				prob_nuc[2][1] = 0.5 - best_error13;
				prob_nuc[2][2] = 0.5 - best_error13;
				prob_nuc[2][3] = best_error13;
				prob_nuc[3][1] = best_error13;
				prob_nuc[3][2] = 1.0 - best_error;
				prob_nuc[3][3] = best_error13;
				// printf("Finished calculating the conditional probabilities\n");
  				for (ig=1;ig<=nsample;ig++){
					// printf("coverage of individual %d: %d\n", ig, ind_coverage[ig]);
					if (ind_coverage[ig] >= min_cov && ind_coverage[ig] <= max_cov) {
						// printf("Passed the coverage requirement\n");
						maxll = -10000000000.0;
						sub_maxll = -10000000000.0;
						id_best_geno = 0;
						ml_ind_read[ig][1] = read_count[ig][n1];
						ml_ind_read[ig][2] = read_count[ig][n2];
						ml_ind_read[ig][3] = ind_coverage[ig] - ml_ind_read[ig][1] - ml_ind_read[ig][2];
						// Calculate the sum of the weights over the genotypes
						sum_weights = 0.0;
						for (gg=1;gg<=3;gg++){
							sum_weights = sum_weights + best_geno_freq[gg]*pow(prob_nuc[gg][1], (double)ml_ind_read[ig][1])*pow(prob_nuc[gg][2], (double)ml_ind_read[ig][2])*pow(prob_nuc[gg][3], (double)ml_ind_read[ig][3]);
						}
						// printf("read counts of individual %d: %d\t%d\t%d\n", ig, ml_ind_read[ig][1], ml_ind_read[ig][2], ml_ind_read[ig][3]);
  						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							// Calculate the probability of the observed set of the site-specific nucleotide reads
							prob_data[ig][ml_geno] = ( best_geno_freq[ml_geno]*pow(prob_nuc[ml_geno][1], (double)ml_ind_read[ig][1])*pow(prob_nuc[ml_geno][2], (double)ml_ind_read[ig][2])*pow(prob_nuc[ml_geno][3], (double)ml_ind_read[ig][3]) )/sum_weights;
							// Take the natural log of the probability
							if (prob_data[ig][ml_geno] > 0.0) {
								llg[ml_geno] = log(prob_data[ig][ml_geno]);
							} else {
								llg[ml_geno] = -10000000001.0;
							}
							// printf("Finished calculating the log-likelihood\n");
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								sub_maxll = maxll;
								maxll = llg[ml_geno];
								if (ml_geno == 1) {	// homozygous for the major allele
									if (ml_ind_read[ig][1] > 0) {
										id_best_geno = 1;
										if (n1 == 1) {		// major allele is A
											best_geno[ig] = "AA";
										} else if (n1 == 2) {		// major allele is C
											best_geno[ig] = "CC";
										} else if (n1 == 3) {          // major allele is G
                                                					best_geno[ig] = "GG";
										} else if (n1 == 4) {          // major allele is T
                                                					best_geno[ig] = "TT";
										}
									}
								} else if (ml_geno == 2) {	// heterozygous
									if (ml_ind_read[ig][1]*ml_ind_read[ig][2] > 0) {
										id_best_geno = 2;
										if (n1*n2 == 2) {	// Alleles are A and C
											best_geno[ig] = "AC";
										} else if (n1*n2 == 3) {     // Alleles are A and G
											best_geno[ig] = "AG";
										} else if (n1*n2 == 4) {     // Alleles are A and T
                                                					best_geno[ig] = "AT";
										} else if (n1*n2 == 6) {     // Alleles are C and G
                                                					best_geno[ig] = "CG";
										} else if (n1*n2 == 8) {     // Alleles are C and T
                                                					best_geno[ig] = "CT";
										} else if (n1*n2 == 12) {     // Alleles are G and T
                                                					best_geno[ig] = "GT";
										}
									}
								} else if (ml_geno == 3) {	// homozygous for the minor allele
									if (ml_ind_read[ig][2] > 0) {
										id_best_geno = 3;
										if (n2 == 1) {         // minor allele is A
                                                					best_geno[ig] = "AA";
                                        					} else if (n2 == 2) {          // minor allele is C
                                                					best_geno[ig] = "CC";
                                        					} else if (n2 == 3) {          // minor allele is G
                                                					best_geno[ig] = "GG";
                                        					} else if (n2 == 4) {          // minor allele is T
                                                					best_geno[ig] = "TT";
                                        					}
									}
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						/*
						if (site == 44) {
							printf("ig: %d\tprob_data[1]: %f\tprob_data[2]: %f\tprob_data[3]: %f\n", ig, prob_data[ig][1], prob_data[ig][2], prob_data[ig][3]);
							printf("major allele: %d\tminor allele: %d\n", n1, n2);
						}
						*/
						if (maxll > sub_maxll) {
							if (id_best_geno == 1) {
								tot_count_geno = tot_count_geno + 1;
								count_geno[1] = count_geno[1] + 1;
								major_allele_count = major_allele_count + 2;
							} else if (id_best_geno == 2) {
								tot_count_geno = tot_count_geno + 1;
								count_geno[2] = count_geno[2] + 1;
								major_allele_count = major_allele_count + 1;
								minor_allele_count = minor_allele_count + 1;
							} else if (id_best_geno == 3) {
								tot_count_geno = tot_count_geno + 1;
								count_geno[3] = count_geno[3] + 1;
								minor_allele_count = minor_allele_count + 2;
							}
							// printf("%d\t%d\t%d\t%d\t%s\t%s\t%f\t%f\t%f\t%f\t%s\t%f\n", ig, ml_ind_read[ig][1], ml_ind_read[ig][2], ml_ind_read[ig][3], major_allele.c_str(), minor_allele.c_str(), llg[1], llg[2], llg[3], max_ll, best_geno.c_str(), best_error);
						} else {
							best_geno[ig] = "NA";
							// printf("%d\t%d\t%d\t%d\t%s\t%s\t%f\t%f\t%f\t%f\tNA\t%f\n", ig, ml_ind_read[ig][1], ml_ind_read[ig][2], ml_ind_read[ig][3], major_allele.c_str(), minor_allele.c_str(), llg[1], llg[2], llg[3], max_ll, best_error);
						}
						if (best_geno[ig] != "NA") {	// genotype call available
							// Filter the genotype call using the binomial test
							if ( anc.substr(0,1) != anc.substr(1,1) ) {	// inferred ancestral genotype heterozygous
								prob_ind_het = 0.0;
								if (ml_ind_read[ig][2] < ml_ind_read[ig][1]) {	// minor nucleotide read less abundant
									for (dg=ml_ind_read[ig][2]; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
								} else if (ml_ind_read[ig][2] == ml_ind_read[ig][1]) {	// nucleotide reads equally abundant
									for (dg=ml_ind_read[ig][2]-1; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
									prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], ml_ind_read[ig][2], (double)0.5-best_error13);
								} else {	// major nucleotide read less abundant
									for (dg=ml_ind_read[ig][1]; dg>=0; dg--) {
										prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
									}
									prob_ind_het = (double)2.0*prob_ind_het;
								}
								if (prob_ind_het < cv_het) {	// unlikely to be heterozygous
									if (ml_ind_read[ig][1] == 0) {
										best_geno[ig] = minor_allele + minor_allele;
										f_minor_allele_count = f_minor_allele_count + 2;
										num_mut_lines = num_mut_lines + 1;
									} else if (ml_ind_read[ig][2] == 0) {
										best_geno[ig] = major_allele + major_allele;
                                                                                f_major_allele_count = f_major_allele_count + 2;
                                                                                num_mut_lines = num_mut_lines + 1;
									} else {
										best_geno[ig] = anc;
										f_count_geno[2] = f_count_geno[2] + 1;
										f_major_allele_count = f_major_allele_count + 1;
										f_minor_allele_count = f_minor_allele_count + 1;
									}
								} else {	// likely to be heterozygous
									best_geno[ig] = anc;
                                                                        f_count_geno[2] = f_count_geno[2] + 1;
                                                                        f_major_allele_count = f_major_allele_count + 1;
                                                                        f_minor_allele_count = f_minor_allele_count + 1;
								}
							} else {	// inferred ancestral genotype homozygous
								base_prob_ind_hom = binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], ml_ind_read[ig][1], (double)1.0-best_error13);
								prob_ind_hom = base_prob_ind_hom;
								for (dg=ml_ind_read[ig][1]+ml_ind_read[ig][2]; dg>=0; dg--) {
									if (binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)1.0-best_error13) < base_prob_ind_hom) {
										prob_ind_hom = prob_ind_hom + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)1.0-best_error13);
									}
								}
								if (prob_ind_hom < cv_hom) {    // unlikely to be major homozygous
									prob_ind_het = 0.0;
									if (ml_ind_read[ig][2] < ml_ind_read[ig][1]) {  // minor nucleotide read less abundant
                                                                        	for (dg=ml_ind_read[ig][2]; dg>=0; dg--) {
                                                                                	prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
                                                                        	}
										prob_ind_het = (double)2.0*prob_ind_het;
									} else if (ml_ind_read[ig][2] == ml_ind_read[ig][1]) {  // nucleotide reads equally abundant
										for (dg=ml_ind_read[ig][2]-1; dg>=0; dg--) {
											prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
										}
										prob_ind_het = (double)2.0*prob_ind_het;
										prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], ml_ind_read[ig][2], (double)0.5-best_error13);
                                                               		} else {
                                                                        	for (dg=ml_ind_read[ig][1]; dg>=0; dg--) {
                                                                                	prob_ind_het = prob_ind_het + binomial(ml_ind_read[ig][1]+ml_ind_read[ig][2], dg, (double)0.5-best_error13);
                                                                        	}
										prob_ind_het = (double)2.0*prob_ind_het;
									}
									if (prob_ind_het > cv_het && ml_ind_read[ig][1] >= 2 && ml_ind_read[ig][2] >= 2) {	// likely to be heterozygous
										best_geno[ig] = major_allele + minor_allele;
										f_major_allele_count = f_major_allele_count + 1;
										f_minor_allele_count = f_minor_allele_count + 1;
										num_mut_lines = num_mut_lines + 1;
									} else {
										if (ml_ind_read[ig][1] == 0) {
											best_geno[ig] = minor_allele + minor_allele;
											f_minor_allele_count = f_minor_allele_count + 2;
											num_mut_lines = num_mut_lines + 1;
										} else {
											best_geno[ig] = anc;
											f_count_geno[1] = f_count_geno[1] + 1;
											f_major_allele_count = f_major_allele_count + 2;
										}
									}
								} else {        // likely to be major homozygous
									best_geno[ig] = anc;
									f_count_geno[1] = f_count_geno[1] + 1;
                                                                        f_major_allele_count = f_major_allele_count + 2;
                                                                }
							}
						}
					}
  				}
  				// printf("scaffold: %s\tposition: %d\tsample size: %d\tallele 1 homo: %f\thetero: %f\tallele 2 homo: %f\n", scaffold.c_str(), site, tot_count_geno, (double)count_geno[1]/tot_count_geno, (double)count_geno[2]/tot_count_geno, (double)count_geno[3]/tot_count_geno);
				// Re-infer the ancestral genotype removing MA lines with mutations
				if (f_count_geno[1] > f_count_geno[2]) {
					f_anc = major_allele + major_allele;
				} else {
					f_anc = major_allele + minor_allele;
				}
  				if (count_geno[1] < tot_count_geno) {
					if ( as == 1 || (f_anc == anc && f_minor_allele_count > 0 && num_mut_lines == 1 ) ) {	// Keep the SNP site only if the ancestral genotype keeps appearing consistent, the site remains polymorphic, and only one MA line has mutations
						if (rn == 0) {
							fprintf(outstream, "%s\t%d\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
							// printf("%s\t%d\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
						} else if (rn == 1) {
							fprintf(outstream, "%s\t%d\t%s\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
							// printf("%s\t%d\t%s\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
						}
						if (gp == 0) {
							for(ig=1;ig<nsample;ig++) {
								fprintf(outstream, "%s\t", best_geno[ig].c_str());
								// printf("%s\t", best_geno[ig].c_str());
							}
							fprintf(outstream, "%s\n", best_geno[ig].c_str());
							// printf("%s\n", best_geno[ig].c_str());
						} else if (gp == 1) {
							for(ig=1;ig<nsample;ig++) {
								if (best_geno[ig] != "NA") {
									fprintf(outstream, "%s\t%f/%f\t", best_geno[ig].c_str(), prob_data[ig][1], prob_data[ig][3]);
									// printf("%s\t%f/%f\t", best_geno[ig].c_str(), prob_data[ig][1], prob_data[ig][3]);
								} else {
									fprintf(outstream, "%s\tNA\t", best_geno[ig].c_str());
									// printf("%s\tNA\t", best_geno[ig].c_str());
								}
							}
							if (best_geno[ig] != "NA") {
								fprintf(outstream, "%s\t%f/%f\n", best_geno[ig].c_str(), prob_data[nsample][1], prob_data[nsample][3]);
								// printf("%s\t%f/%f\n", best_geno[ig].c_str(), prob_data[nsample][1], prob_data[nsample][3]);
							} else {
								fprintf(outstream, "%s\tNA\n", best_geno[ig].c_str());
								// printf("%s\tNA\n", best_geno[ig].c_str());
							}
						}
  					}
				} else {
					if (as == 1) {	// Print out the genotype calls at monomorphic sites
						if (rn == 0) {
							fprintf(outstream, "%s\t%d\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
							// printf("%s\t%d\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
						} else if (rn == 1) {
							fprintf(outstream, "%s\t%d\t%s\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
							// printf("%s\t%d\t%s\t%f\t%d\t%s\t%d\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage, f_anc.c_str(), tot_count_geno, f_minor_allele_count);
						}
						if (gp == 0) {
                                                        for(ig=1;ig<nsample;ig++) {
                                                                fprintf(outstream, "%s\t", best_geno[ig].c_str());
                                                                // printf("%s\t", best_geno[ig].c_str());
                                                        }
                                                        fprintf(outstream, "%s\n", best_geno[ig].c_str());
                                                        // printf("%s\n", best_geno[ig].c_str());
                                                } else if (gp == 1) {
                                                        for(ig=1;ig<nsample;ig++) {
                                                                if (best_geno[ig] != "NA") {
                                                                        fprintf(outstream, "%s\t%f/%f\t", best_geno[ig].c_str(), prob_data[ig][1], prob_data[ig][3]);
                                                                        // printf("%s\t%f/%f\t", best_geno[ig].c_str(), prob_data[ig][1], prob_data[ig][3]);
                                                                } else {
                                                                        fprintf(outstream, "%s\tNA\t", best_geno[ig].c_str());
                                                                        // printf("%s\tNA\t", best_geno[ig].c_str());
                                                                }
                                                        }
                                                        if (best_geno[ig] != "NA") {
                                                                fprintf(outstream, "%s\t%f/%f\n", best_geno[ig].c_str(), prob_data[nsample][1], prob_data[nsample][3]);
                                                                // printf("%s\t%f/%f\n", best_geno[ig].c_str(), prob_data[nsample][1], prob_data[nsample][3]);
                                                        } else {
                                                                fprintf(outstream, "%s\tNA\n", best_geno[ig].c_str());
                                                                // printf("%s\tNA\n", best_geno[ig].c_str());
                                                        }
                                                }
					}
				}
			} else {	// Print out the genotype calls at monomorphic sites
				if (rn == 0) {
					fprintf(outstream, "%s\t%d\t%f\t%d\t", scaffold.c_str(), site, best_error, pop_coverage);
                			// printf("%s\t%d\t%f\t%d\t", scaffold.c_str(), site, best_error, pop_coverage);
				} else if (rn == 1) {
					fprintf(outstream, "%s\t%d\t%s\t%f\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage);
					// printf("%s\t%d\t%s\t%f\t%d\t", scaffold.c_str(), site, ref_nuc.c_str(), best_error, pop_coverage);
				}
				// The ancestral genotype inferred to be the major homozygote
				if (n1 == 1) {          // major allele is A
					anc = "AA";
				} else if (n1 == 2) {           // major allele is C
					anc = "CC";
				} else if (n1 == 3) {          // major allele is G
					anc = "GG";
				} else if (n1 == 4) {          // major allele is T
					anc = "TT";
				}
				tot_count_geno = 0;
				f_minor_allele_count = 0;
				for (ig=1;ig<=nsample;ig++){
					// printf("coverage of individual %d: %d\n", ig, ind_coverage[ig]);
                    			if (ind_coverage[ig] >= min_cov && ind_coverage[ig] <= max_cov) {
						tot_count_geno = tot_count_geno + 1;
						if (n1 == 1) {          // major allele is A
							best_geno[ig] = "AA";
						} else if (n1 == 2) {           // major allele is C
							best_geno[ig] = "CC";
						} else if (n1 == 3) {          // major allele is G
							best_geno[ig] = "GG";
						} else if (n1 == 4) {          // major allele is T
							best_geno[ig] = "TT";
						}
					} else {
						best_geno[ig] = "NA";
					}
				}
				fprintf(outstream, "%s\t%d\t%d\t", anc.c_str(), tot_count_geno, f_minor_allele_count);
				// printf("%s\t%d\t%d\t", anc.c_str(), tot_count_geno, f_minor_allele_count);
				for (ig=1;ig<=nsample;ig++){
					if (gp == 0) {
						if (ig < nsample) {
							fprintf(outstream, "%s\t", best_geno[ig].c_str());
							// printf("%s\t", best_geno[ig].c_str());
						} else {
							fprintf(outstream, "%s\n", best_geno[ig].c_str());
							// printf(outstream, "%s\n", best_geno[ig].c_str());
						}
					} else if (gp == 1) {
						if (ig < nsample) {
							fprintf(outstream, "%s\tNA\t", best_geno[ig].c_str());
							// printf("%s\tNA\t", best_geno[ig].c_str());
						} else {
							fprintf(outstream, "%s\tNA\n", best_geno[ig].c_str());
							// printf(outstream, "%s\tNA\n", best_geno[ig].c_str());
						}
					}
				}
			}
		} else {
			if (as == 1) {	// Print out the sites without ML estimates
				if (rn == 0) {
					fprintf(outstream, "%s\t%d\tNA\t%d\tNA\tNA\tNA\t", scaffold.c_str(), site, pop_coverage);
					// printf("%s\t%d\tNA\t%d\tNA\tNA\tNA\t", scaffold.c_str(), site, pop_coverage);
				} else if (rn == 1) {
					fprintf(outstream, "%s\t%d\t%s\tNA\t%d\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
					// printf("%s\t%d\t%s\tNA\t%d\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
				}
				if (gp == 0) {
					for (ig=1;ig<=nsample;ig++){
						if (ig < nsample) {
							fprintf(outstream, "NA\t");
                    					// printf(outstream, "NA\t");
						} else {
							fprintf(outstream, "NA\n");
                    					// printf(outstream, "NA\n");
                				}
					}
            			} else if (gp == 1) {
					for (ig=1;ig<=nsample;ig++){
						if (ig < nsample) {
							fprintf(outstream, "NA\tNA\t");
                    					// printf(outstream, "NA\tNA\t");
                				} else {
							fprintf(outstream, "NA\tNA\n");
                    					// printf(outstream, "NA\tNA\n");
                				}
            				}
				}
			}
		}
  	}

	// Close the input and output files
	inputFile.close();
	fclose(outstream);

  	return (0);
}
