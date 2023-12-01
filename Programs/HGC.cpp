/* Updated on 03/29/16 */

/*
 
Program HGC.cpp to call the genotype of each individual from pro files 
of individual high-throughput sequencing data of multiple diploid individuals 
from a population.  In this estimator, the assumption of at most two alleles 
segregating at a site is relaxed.  In addition, called genotypes are reported 
only at significantly polymorphic sites in this version.  Also, the number of 
alleles are based on significant variants according to the likelihood-ratio test.
The statistical significance of alleles in heterozygotes is examined more conservatively 
in this version.  The mean error rate among individuals is reported at each site in 
this version.  The allele counts are reported in this version.  The maximum allowed 
error rate can be specified by users in this version.  In this version, information 
at every site of the reference sequence is shown in the output.  

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

int main(int argc, char *argv[]){

	// Default values of the options
	const char* in_file_name = {"In_HGC.txt"};
	const char* out_file_name = {"Out_HGC.txt"};
	int min_cov = 1;	// minimum coverage cut-off required to call a genotype
	int max_cov = 2000000000;	// maximum coverage cut-off allowed to call a genotype
	double cv = 3.841;
	double max_e = 0.5;
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
		} else if (strcmp(argv[argz], "-cv") == 0) {
			sscanf(argv[++argz], "%lf", &cv);
		} else if (strcmp(argv[argz], "-max_e") == 0) {
			sscanf(argv[++argz], "%lf", &max_e);
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
		fprintf(stderr, "       -cv <f>: specify the chi-square critical value for the polymorphism test\n");
		fprintf(stderr, "       -max_e <f>: specify the maximum error-rate estimate allowed\n");
		exit(1);
	}
	
	if (min_cov > 1) {
		printf("min_cov: %d\n", min_cov);
	}
	if (max_cov < 2000000000) {
		printf("max_cov: %d\n", max_cov);
	}
	if (max_e != 0.5) {
		printf("max_e: %f\n", max_e);
	}

	ifstream inputFile(in_file_name);	// Try to open the input file
	if ( !inputFile.is_open() ) {	// Exit on failure
		printf("Cannot open %s for reading.\n", in_file_name);
		exit(1);
	}

	// Read the header
	string line;    // String buffer
	string h_scaf, h_site, h_ref_nuc;	// Stores header names	
	getline(inputFile, line);
	istringstream ss(line);
        ss >> h_scaf >> h_site >> h_ref_nuc;
	string str;	// Temporarily stores each individual ID
	vector <string> id_ind;      // Stores individual IDs
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
	
  	int ag, ig, gg, jg, kg, lg, mg, ng, size_quartet, digit, count_nuc;
	string scaffold;
  	int site;
  	string ref_nuc;
  	string quartet[nsample+1];
  	int ind_coverage[nsample+1], num_digit[5], read_count[nsample+1][5], pop_read[5];
	int pop_coverage, n1, n2;
	string major_allele, minor_allele, best_geno[nsample+1];
	string test_geno;	// temporarily genotype for examining the stasistical significance of the alleles
	int count_geno[11], tot_count_geno, count_allele[5], ml_geno;
	double maxll, sub_maxll, error, llg[11], best_error;
	double null_error, null_maxll, ind_llstat, max_ind_llstat[5];
	double sum_best_error, mean_best_error;
	int num_alleles;
	int a1, a2, a3, a4; 

	// point to the output file
	FILE *outstream;
	
  	// Open the output file 
  	outstream = fopen( out_file_name, "w" );

  	// Print the field names
	fprintf(outstream, "scaffold\tposition\tref_nuc\tnum_alleles\tmajor_allele\tminor_allele\tpop_coverage\tsample_size\tmajor_allele_count\tminor_allele_count\tmean_best_error\t");
	// printf("scaffold\tposition\tref_nuc\tnum_alleles\tmajor_allele\tminor_allele\tpop_coverage\tsample_size\tmajor_allele_count\tminor_allele_count\tmean_best_error\t");
	for(ig=0;ig<nsample-1;ig++) {
		fprintf(outstream, "%s\t", id_ind[ig].c_str());
        	// printf("%s\t", id_ind[ig].c_str());
	}
	fprintf(outstream, "%s\tmax_ind_llstat[1]\tmax_ind_llstat[2]\tmax_ind_llstat[3]\tmax_ind_llstat[4]\tcount_allele[a1]\tcount_allele[a2]\tcount_allele[a3]\tcount_allele[a4]\n", id_ind[nsample-1].c_str());
	// printf("%s\tmax_ind_llstat[1]\tmax_ind_llstat[2]\tmax_ind_llstat[3]\tmax_ind_llstat[4]\tcount_allele[a1]\tcount_allele[a2]\tcount_allele[a3]\tcount_allele[a4]\n", id_ind[nsample-1].c_str());
  	// Read the main data 
  	while ( getline(inputFile, line) ) {
  		istringstream ss(line);
		ss >> scaffold >> site >> ref_nuc;
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
			printf("%s\t%d\t%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str(), ig, id_ind[ig-1].c_str());
			fprintf(outstream, "%s\t%d\t%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str(), ig, id_ind[ig-1].c_str());
			for (jg = 1; jg <= 4; jg++) {
				printf("%d\t", read_count[ig][jg]);
				fprintf(outstream, "%d\t", read_count[ig][jg]);
			}
			printf("%d\n", ind_coverage[ig]);
			fprintf(outstream, "%d\n", ind_coverage[ig]);
		}
		*/
		// ML estimation starts here
        	// Calculate the number of four different nucleotide reads in the population sample.
        	pop_coverage = 0;
        	for (jg = 1; jg <= 4; jg++) {
                	pop_read[jg] = 0;
                	for (ig = 1; ig <= nsample; ig++) {
                        	pop_read[jg] = pop_read[jg] + read_count[ig][jg];
                	}
                	pop_coverage = pop_coverage + pop_read[jg];
        	}
        	// printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", scaffold.c_str(), site, pop_read[1], pop_read[2], pop_read[3], pop_read[4], pop_coverage);
        	// fprintf(outstream, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", scaffold.c_str(), site, pop_read[1], pop_read[2], pop_read[3], pop_read[4], pop_coverage);
        	// Find the two most abundant nucleotides at the site
        	if (pop_read[1] >= pop_read[2]) {
                	n1 = 1;
                	n2 = 2;
        	} else {
                	n1 = 2;
                	n2 = 1;
        	}
        	if (pop_read[3] > pop_read[n1]) {
                	n2 = n1;
                	n1 = 3;
        	} else if (pop_read[3] > pop_read[n2]) {
			n2 = 3;
		}
		if (pop_read[4] > pop_read[n1]) {
                	n2 = n1;
                	n1 = 4;
		} else if (pop_read[4] > pop_read[n2]) {
			n2 = 4;
		}
		// printf("%s\t%d\t%d\t%d\n", scaffold.c_str(), site, n1, n2);
		// fprintf(outstream, "%s\t%d\t%d\t%d\n", scaffold.c_str(), site, n1, n2);
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
		// The default individual genotypes are NA
		for (ig=1;ig<=nsample;ig++) {
			best_geno[ig] = "NA";
			// printf("%s\n", best_geno[ig].c_str());
		}
		// printf("Just before the ML estimation\n");
		// Carry out the ML estimation only when there are at least two different nucleotides in the population sample	 
		if (pop_read[n1] < pop_coverage) {
			tot_count_geno = 0;
			for (gg=1;gg<=10;gg++) {
				count_geno[gg] = 0;
			}
			for (ag=1;ag<=4;ag++) {
				count_allele[ag] = 0;
				max_ind_llstat[ag] = 0.0;
			}
			sum_best_error = 0.0;
  			for (ig=1;ig<=nsample;ig++){
				// printf("coverage of individual %d: %d\n", ig, ind_coverage[ig]); 
				if (ind_coverage[ig] >= min_cov && ind_coverage[ig] <= max_cov) {
					// printf("Passed the coverage requirement\n");
					maxll = -10000000000.0;
					sub_maxll = -10000000000.0;
					// printf("read counts of individual %d: %d\t%d\t%d\t%d\n", ig, read_count[ig][1], read_count[ig][2], read_count[ig][3], read_count[ig][4]);
					if (read_count[ig][1] == ind_coverage[ig]) {
						best_error = 0.0;
						best_geno[ig] = "AA";
						maxll = 0.0;
						tot_count_geno = tot_count_geno + 1;
						count_geno[1] = count_geno[1] + 1;
						count_allele[1] = count_allele[1] + 2;
						null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
						if (n1 == 1) {
							null_maxll = 0.0;
						} else {
							null_maxll = read_count[ig][1]*log(null_error/3.0);
						}
						ind_llstat = 2.0*(maxll-null_maxll);
						if (ind_llstat > max_ind_llstat[1]) {
							max_ind_llstat[1] = ind_llstat;
						}
					} else if (read_count[ig][2] == ind_coverage[ig]) {
						best_error = 0.0;
						best_geno[ig] = "CC";
						maxll = 0.0;
						tot_count_geno = tot_count_geno + 1;
						count_geno[5] = count_geno[5] + 1;
						count_allele[2] = count_allele[2] + 2;
						null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
						if (n1 == 2) {
							null_maxll = 0.0;
						} else {
							null_maxll = read_count[ig][2]*log(null_error/3.0);
						}
						ind_llstat = 2.0*(maxll-null_maxll);
						if (ind_llstat > max_ind_llstat[2]) {
							max_ind_llstat[2] = ind_llstat;
						}
					} else if (read_count[ig][3] == ind_coverage[ig]) {
						best_error = 0.0;
						best_geno[ig] = "GG";
						maxll = 0.0;
						tot_count_geno = tot_count_geno + 1;
						count_geno[8] = count_geno[8] + 1;
						count_allele[3] = count_allele[3] + 2;
						null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
						if (n1 == 3) {
							null_maxll = 0.0;
						} else {
							null_maxll = read_count[ig][3]*log(null_error/3.0);
						}
						ind_llstat = 2.0*(maxll-null_maxll);
						if (ind_llstat > max_ind_llstat[3]) {
							max_ind_llstat[3] = ind_llstat;
						}
					} else if (read_count[ig][4] == ind_coverage[ig]) {
						best_error = 0.0;
						best_geno[ig] = "TT";
						maxll = 0.0;
						tot_count_geno = tot_count_geno + 1;
						count_geno[10] = count_geno[10] + 1;
						count_allele[4] = count_allele[4] + 2;
						null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
						if (n1 == 4) {
							null_maxll = 0.0;
						} else {
							null_maxll = read_count[ig][4]*log(null_error/3.0);
						}
						ind_llstat = 2.0*(maxll-null_maxll);
						if (ind_llstat > max_ind_llstat[4]) {
							max_ind_llstat[4] = ind_llstat;
						}
					} else if ( (read_count[ig][1]+read_count[ig][2]) == ind_coverage[ig] ) { // nucleotide reads A and C only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype AA
								error = read_count[ig][2]/(double)ind_coverage[ig];
								if (error > max_e) {
									error = max_e;
								}
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + read_count[ig][2]*log(error/3.0); 		
							} else if (ml_geno == 2){	// candidate genotype AC
								error = 0.0;
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][2])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype CC
								error = read_count[ig][1]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + read_count[ig][1]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AC";
								} else if (ml_geno == 3) {	// CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "CC";
								}			 
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}	// end of the loop over the candidate genotypes
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AC") {
								count_geno[2] = count_geno[2] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[2] = count_allele[2] + 1;
							} else if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + read_count[ig][2]*log(null_error/3.0);
							} else if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + read_count[ig][1]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][2])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AC") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][2]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AC";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "AC";
										}
									}
								} else {
									test_geno = "AC";
								}
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][3]) == ind_coverage[ig] ) { // nucleotide reads A and G only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							if (ml_geno == 1) {	// candidate genotype AA
								error = read_count[ig][3]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + read_count[ig][3]*log(error/3.0);
							} else if (ml_geno == 2){	// candidate genotype AG
								error = 0.0;
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][3])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype GG
								error = read_count[ig][1]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + read_count[ig][1]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AG";
								} else if (ml_geno == 3) {	// GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "GG";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AG") {
								count_geno[3] = count_geno[3] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + read_count[ig][3]*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + read_count[ig][1]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][3])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AG") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][3]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AG";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "GG";
										} else {
											test_geno = "AG";
										}
									}
								} else {
									test_geno = "AG";
								}
							} else {
								test_geno = best_geno[ig];
							}	
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads A and T only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							if (ml_geno == 1) {	// candidate genotype AA
								error = read_count[ig][4]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 2){	// candidate genotype AT
								error = 0.0;
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][4])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype TT
								error = read_count[ig][1]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + read_count[ig][1]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "AT";
								} else if (ml_geno == 3) {	// TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AT") {
								count_geno[4] = count_geno[4] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + read_count[ig][4]*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + read_count[ig][1]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AT") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][4]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AT";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "TT";
										} else {
											test_geno = "AT";
										}
									}
								} else {
									test_geno = "AT";
								}
							} else {
								test_geno = best_geno[ig];
							}
							/*
							if (ig == 53) {
								printf("test_geno: %s\tind_llstat: %f\n", test_geno.c_str(), ind_llstat);
							}
							*/			
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][2]+read_count[ig][3]) == ind_coverage[ig] ) { // nucleotide reads C and G only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							if (ml_geno == 1) {	// candidate genotype CC
								error = read_count[ig][3]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + read_count[ig][3]*log(error/3.0);
								/*
								if (ig == 26) {
									printf("n_C: %d\tn_G: %d\terror: %f\tLL(CC): %f\n", read_count[ig][2], read_count[ig][3], error, llg[ml_geno]);
								}
								*/
							} else if (ml_geno == 2){	// candidate genotype CG
								error = 0.0;
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][3])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype GG
								error = read_count[ig][2]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + read_count[ig][2]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "CC";
								} else if (ml_geno == 2) {	// CG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "CG";
								} else if (ml_geno == 3) {	// GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "GG";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CG") {
								count_geno[6] = count_geno[6] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + read_count[ig][3]*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + read_count[ig][2]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][2]+read_count[ig][3])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "CG") {
								if (ind_llstat > cv) {
									if (read_count[ig][2] >= read_count[ig][3]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "CG";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "GG";
										} else {
											test_geno = "CG";
										}
									}
								} else {
									test_geno = "CG";
								}
							} else {
								test_geno = best_geno[ig];
							}		
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									}
								}
							}	
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][2]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads C and T only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							if (ml_geno == 1) {	// candidate genotype CC
								error = read_count[ig][4]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 2){	// candidate genotype CT
								error = 0.0;
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][4])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype TT
								error = read_count[ig][2]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + read_count[ig][2]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "CC";
								} else if (ml_geno == 2) {	// CT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "CT";
								} else if (ml_geno == 3) {	// TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CT") {
								count_geno[7] = count_geno[7] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + read_count[ig][4]*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + read_count[ig][2]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][2]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "CT") {
								if (ind_llstat > cv) {
									if (read_count[ig][2] >= read_count[ig][4]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "CT";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "TT";
										} else {
											test_geno = "CT";
										}
									}
								} else {
									test_geno = "CT";
								}
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][3]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads G and T only
						// Loop over three candidate genotypes of the individual
						for (gg=1;gg<=3;gg++){
							ml_geno = gg;
							if (ml_geno == 1) {	// candidate genotype GG
								error = read_count[ig][4]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 2){	// candidate genotype GT
								error = 0.0;
								llg[ml_geno] = (read_count[ig][3]+read_count[ig][4])*log(0.5-error/3.0);
							} else if (ml_geno == 3){	// candidate genotype TT
								error = read_count[ig][3]/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + read_count[ig][3]*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "GG";
								} else if (ml_geno == 2) {	// GT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "GT";
								} else if (ml_geno == 3) {	// TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error;
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							} else if (best_geno[ig] == "GT") {
								count_geno[9] = count_geno[9] + 1;
								count_allele[3] = count_allele[3] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + read_count[ig][4]*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + read_count[ig][3]*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][3]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "GT") {
								if (ind_llstat > cv) {
									if (read_count[ig][3] >= read_count[ig][4]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "GG";
										} else {
											test_geno = "GT";
										}
									} else {
										if (maxll - llg[3] <= cv) {
											test_geno = "TT";
										} else {
											test_geno = "GT";
										}
									}
								} else {
									test_geno = "GT";
								}
							} else {
								test_geno = best_geno[ig];
							}		
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][2]+read_count[ig][3]) == ind_coverage[ig] ) { // nucleotide reads A, C and G only
							// Loop over six candidate genotypes of the individual
						for (gg=1;gg<=6;gg++) {
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype AA
								error = (ind_coverage[ig]-read_count[ig][1])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + (read_count[ig][2]+read_count[ig][3])*log(error/3.0); 		
							} else if (ml_geno == 2) {	// candidate genotype AC
								if ( 1.5*( read_count[ig][3]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][3]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][2])*log(0.5-error/3.0) + read_count[ig][3]*log(error/3.0);
							} else if (ml_geno == 3) {	// candidate genotype AG
								if ( 1.5*( read_count[ig][2]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][2]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][3])*log(0.5-error/3.0) + read_count[ig][2]*log(error/3.0);
							} else if (ml_geno == 4) {	// Candidate genotype CC
								error = (ind_coverage[ig]-read_count[ig][2])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + (read_count[ig][1]+read_count[ig][3])*log(error/3.0);
							} else if (ml_geno == 5) {	// Candidate genotype CG
								if ( 1.5*( read_count[ig][1]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][1]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][3])*log(0.5-error/3.0) + read_count[ig][1]*log(error/3.0);
							} else if (ml_geno == 6) {	// Candidate genotype GG
								error = (ind_coverage[ig]-read_count[ig][3])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + (read_count[ig][1]+read_count[ig][2])*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AC 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AC";
								} else if (ml_geno == 3) {	// AG 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AG";
								} else if (ml_geno == 4) { // CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CC";
								} else if (ml_geno == 5) { // CG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CG";
								} else if (ml_geno == 6) { // GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GG";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {  
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AC") {
								count_geno[2] = count_geno[2] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[2] = count_allele[2] + 1;
							} else if (best_geno[ig] == "AG") {
								count_geno[3] = count_geno[3] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CG") {
								count_geno[6] = count_geno[6] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + (read_count[ig][2]+read_count[ig][3])*log(null_error/3.0);
							} else if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][3])*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][2])*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][2]+read_count[ig][3])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AC") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][2]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AC";
										}
									} else {
										if (maxll - llg[4] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "AC";
										}
									}
								} else {
									test_geno = "AC";
								}
							} else if (best_geno[ig] == "AG") {
								if (ind_llstat > cv) {
                                                                        if (read_count[ig][1] >= read_count[ig][3]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "AA";
                                                                                } else {
                                                                                        test_geno = "AG";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "AG";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "AG";
                                                                }
							} else if (best_geno[ig] == "CG") {
								if (ind_llstat > cv) {
                                                                        if (read_count[ig][2] >= read_count[ig][3]) {
                                                                                if (maxll - llg[4] <= cv) {
                                                                                        test_geno = "CC";
                                                                                } else {
                                                                                        test_geno = "CG";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "CG";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "CG";
                                                                }
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][2]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads A, C, and T only
							// Loop over six candidate genotypes of the individual
						for (gg=1;gg<=6;gg++) {
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype AA
								error = (ind_coverage[ig]-read_count[ig][1])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + (read_count[ig][2]+read_count[ig][4])*log(error/3.0); 		
							} else if (ml_geno == 2) {	// candidate genotype AC
								if ( 1.5*( read_count[ig][4]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][4]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][2])*log(0.5-error/3.0) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 3) {	// candidate genotype AT
								if ( 1.5*( read_count[ig][2]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][2]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][2]*log(error/3.0);
							} else if (ml_geno == 4) {	// Candidate genotype CC
								error = (ind_coverage[ig]-read_count[ig][2])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + (read_count[ig][1]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 5) {	// Candidate genotype CT
								if ( 1.5*( read_count[ig][1]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][1]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][1]*log(error/3.0);
							} else if (ml_geno == 6) {	// Candidate genotype TT
								error = (ind_coverage[ig]-read_count[ig][4])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + (read_count[ig][1]+read_count[ig][2])*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AC 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AC";
								} else if (ml_geno == 3) {	// AT 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AT";
								} else if (ml_geno == 4) { // CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CC";
								} else if (ml_geno == 5) { // CT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CT";
								} else if (ml_geno == 6) { // TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {  
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AC") {
								count_geno[2] = count_geno[2] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[2] = count_allele[2] + 1;
							} else if (best_geno[ig] == "AT") {
								count_geno[4] = count_geno[4] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CT") {
								count_geno[7] = count_geno[7] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + (read_count[ig][2]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][2])*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][2]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AC") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][2]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AC";
										}
									} else {
										if (maxll - llg[4] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "AC";
										}
									}
								} else {
									test_geno = "AC";
								}
							} else if (best_geno[ig] == "AT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][1] >= read_count[ig][4]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "AA";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "AT";
                                                                }
							} else if (best_geno[ig] == "CT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][2] >= read_count[ig][4]) {
                                                                                if (maxll - llg[4] <= cv) {
                                                                                        test_geno = "CC";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "CT";
                                                                }
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][3]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads A, G, and T only
							// Loop over six candidate genotypes of the individual
						for (gg=1;gg<=6;gg++) {
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype AA
								error = (ind_coverage[ig]-read_count[ig][1])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + (read_count[ig][3]+read_count[ig][4])*log(error/3.0); 		
							} else if (ml_geno == 2) {	// candidate genotype AG
								if ( 1.5*( read_count[ig][4]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][4]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][3])*log(0.5-error/3.0) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 3) {	// candidate genotype AT
								if ( 1.5*( read_count[ig][3]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][3]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][3]*log(error/3.0);
							} else if (ml_geno == 4) {	// Candidate genotype GG
								error = (ind_coverage[ig]-read_count[ig][3])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + (read_count[ig][1]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 5) {	// Candidate genotype GT
								if ( 1.5*( read_count[ig][1]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][1]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][3]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][1]*log(error/3.0);
							} else if (ml_geno == 6) {	// Candidate genotype TT
								error = (ind_coverage[ig]-read_count[ig][4])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + (read_count[ig][1]+read_count[ig][3])*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AG 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AG";
								} else if (ml_geno == 3) {	// AT 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AT";
								} else if (ml_geno == 4) { // GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GG";
								} else if (ml_geno == 5) { // GT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GT";
								} else if (ml_geno == 6) { // TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {  
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AG") {
								count_geno[3] = count_geno[3] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "AT") {
								count_geno[4] = count_geno[4] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							} else if (best_geno[ig] == "GT") {
								count_geno[9] = count_geno[9] + 1;
								count_allele[3] = count_allele[3] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + (read_count[ig][3]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + (read_count[ig][1]+read_count[ig][3])*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][1]+read_count[ig][3]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AG") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][3]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AG";
										}
									} else {
										if (maxll - llg[4] <= cv) {
											test_geno = "GG";
										} else {
											test_geno = "AG";
										}
									}
								} else {
									test_geno = "AG";
								}
							} else if (best_geno[ig] == "AT") {
								if (ind_llstat > cv) {
                                                                        if (read_count[ig][1] >= read_count[ig][4]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "AA";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "AT";
                                                                }
							} else if (best_geno[ig] == "GT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][3] >= read_count[ig][4]) {
                                                                                if (maxll - llg[4] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "GT";
                                                                }
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][2]+read_count[ig][3]+read_count[ig][4]) == ind_coverage[ig] ) { // nucleotide reads C, G, and T only
							// Loop over six candidate genotypes of the individual
						for (gg=1;gg<=6;gg++) {
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype CC
								error = (ind_coverage[ig]-read_count[ig][2])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + (read_count[ig][3]+read_count[ig][4])*log(error/3.0); 		
							} else if (ml_geno == 2) {	// candidate genotype CG
								if ( 1.5*( read_count[ig][4]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][4]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][3])*log(0.5-error/3.0) + read_count[ig][4]*log(error/3.0);
							} else if (ml_geno == 3) {	// candidate genotype CT
								if ( 1.5*( read_count[ig][3]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][3]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][3]*log(error/3.0);
							} else if (ml_geno == 4) {	// Candidate genotype GG
								error = (ind_coverage[ig]-read_count[ig][3])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + (read_count[ig][2]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 5) {	// Candidate genotype GT
								if ( 1.5*( read_count[ig][2]/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( read_count[ig][2]/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][3]+read_count[ig][4])*log(0.5-error/3.0) + read_count[ig][2]*log(error/3.0);
							} else if (ml_geno == 6) {	// Candidate genotype TT
								error = (ind_coverage[ig]-read_count[ig][4])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + (read_count[ig][2]+read_count[ig][3])*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// CC 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CC";
								} else if (ml_geno == 2) {	// CG 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CG";
								} else if (ml_geno == 3) {	// CT 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CT";
								} else if (ml_geno == 4) { // GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GG";
								} else if (ml_geno == 5) { // GT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GT";
								} else if (ml_geno == 6) { // TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {  
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CG") {
								count_geno[6] = count_geno[6] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "CT") {
								count_geno[7] = count_geno[7] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							} else if (best_geno[ig] == "GT") {
								count_geno[9] = count_geno[9] + 1;
								count_allele[3] = count_allele[3] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + (read_count[ig][3]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + (read_count[ig][2]+read_count[ig][4])*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + (read_count[ig][2]+read_count[ig][3])*log(null_error/3.0);
							} else {
								null_maxll = (read_count[ig][2]+read_count[ig][3]+read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "CG") {
								if (ind_llstat > cv) {
									if (read_count[ig][2] >= read_count[ig][3]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "CG";
										}
									} else {
										if (maxll - llg[4] <= cv) {
											test_geno = "GG";
										} else {
											test_geno = "CG";
										}
									}
								} else {
									test_geno = "CG";
								}
							} else if (best_geno[ig] == "CT") {
								if (ind_llstat > cv) {
                                                                        if (read_count[ig][2] >= read_count[ig][4]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "CC";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "CT";
                                                                }
							} else if (best_geno[ig] == "GT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][3] >= read_count[ig][4]) {
                                                                                if (maxll - llg[4] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[6] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "GT";
                                                                }
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					} else if ( (read_count[ig][1]+read_count[ig][2]+read_count[ig][3]+read_count[ig][4]) == ind_coverage[ig] ) { // All four nucleotide reads present
						// Loop over ten candidate genotypes of the individual
						for (gg=1;gg<=10;gg++) {
							ml_geno = gg; 
							if (ml_geno == 1) {	// candidate genotype AA
								error = (ind_coverage[ig]-read_count[ig][1])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][1]*log(1.0-error) + (read_count[ig][2]+read_count[ig][3]+read_count[ig][4])*log(error/3.0); 		
							} else if (ml_geno == 2) {	// candidate genotype AC
								if ( 1.5*( (read_count[ig][3]+read_count[ig][4])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][3]+read_count[ig][4])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][2])*log(0.5-error/3.0) + (read_count[ig][3]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 3) {	// candidate genotype AG
								if ( 1.5*( (read_count[ig][2]+read_count[ig][4])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][2]+read_count[ig][4])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][3])*log(0.5-error/3.0) + (read_count[ig][2]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 4) {	// candidate genotype AT
								if ( 1.5*( (read_count[ig][2]+read_count[ig][3])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][2]+read_count[ig][3])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][1]+read_count[ig][4])*log(0.5-error/3.0) + (read_count[ig][2]+read_count[ig][3])*log(error/3.0);
							} else if (ml_geno == 5) {	// Candidate genotype CC
								error = (ind_coverage[ig]-read_count[ig][2])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][2]*log(1.0-error) + (read_count[ig][1]+read_count[ig][3]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 6) {	// Candidate genotype CG
								if ( 1.5*( (read_count[ig][1]+read_count[ig][4])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][1]+read_count[ig][4])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][3])*log(0.5-error/3.0) + (read_count[ig][1]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 7) {	// Candidate genotype CT
								if ( 1.5*( (read_count[ig][1]+read_count[ig][3])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][1]+read_count[ig][3])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][2]+read_count[ig][4])*log(0.5-error/3.0) + (read_count[ig][1]+read_count[ig][3])*log(error/3.0);
							} else if (ml_geno == 8) {	// Candidate genotype GG
								error = (ind_coverage[ig]-read_count[ig][3])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][3]*log(1.0-error) + (read_count[ig][1]+read_count[ig][2]+read_count[ig][4])*log(error/3.0);
							} else if (ml_geno == 9) {	// Candidate genotype GT
								if ( 1.5*( (read_count[ig][1]+read_count[ig][2])/(double)ind_coverage[ig] ) > 1.0 ) {
									error = 1.0;
								} else {
									error = 1.5*( (read_count[ig][1]+read_count[ig][2])/(double)ind_coverage[ig] );
								}
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = (read_count[ig][3]+read_count[ig][4])*log(0.5-error/3.0) + (read_count[ig][1]+read_count[ig][2])*log(error/3.0); 
							} else if (ml_geno == 10) {	// Candidate genotype TT
								error = (ind_coverage[ig]-read_count[ig][4])/(double)ind_coverage[ig];
								if (error > max_e) {
                                                                        error = max_e;
                                                                }
								llg[ml_geno] = read_count[ig][4]*log(1.0-error) + (read_count[ig][1]+read_count[ig][2]+read_count[ig][3])*log(error/3.0);
							}
							// Examine whether this is a new ML solution for the individual  
							if (llg[ml_geno] > maxll){
								if (ml_geno == 1) {	// AA 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AA";
								} else if (ml_geno == 2) {	// AC 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AC";
								} else if (ml_geno == 3) {	// AG 
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AG";
								} else if (ml_geno == 4) { // AT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "AT";
								} else if (ml_geno == 5) { // CC
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CC";
								} else if (ml_geno == 6) { // CG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CG";
								} else if (ml_geno == 7) { // CT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "CT";
								} else if (ml_geno == 8) { // GG
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GG";
								} else if (ml_geno == 9) { // GT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "GT";
								} else if (ml_geno == 10) { // TT
									sub_maxll = maxll;
									maxll = llg[ml_geno];
									best_error = error; 
									best_geno[ig] = "TT";
								}
							} else if (llg[ml_geno] == maxll) {
								sub_maxll = maxll;
								maxll = llg[ml_geno];
							} else if (llg[ml_geno] > sub_maxll) {
								sub_maxll = llg[ml_geno];
							}
						}
						if (maxll > sub_maxll) {  
							tot_count_geno = tot_count_geno + 1;
							if (best_geno[ig] == "AA") {
								count_geno[1] = count_geno[1] + 1;
								count_allele[1] = count_allele[1] + 2;
							} else if (best_geno[ig] == "AC") {
								count_geno[2] = count_geno[2] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[2] = count_allele[2] + 1;
							} else if (best_geno[ig] == "AG") {
								count_geno[3] = count_geno[3] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "AT") {
								count_geno[4] = count_geno[4] + 1;
								count_allele[1] = count_allele[1] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "CC") {
								count_geno[5] = count_geno[5] + 1;
								count_allele[2] = count_allele[2] + 2;
							} else if (best_geno[ig] == "CG") {
								count_geno[6] = count_geno[6] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[3] = count_allele[3] + 1;
							} else if (best_geno[ig] == "CT") {
								count_geno[7] = count_geno[7] + 1;
								count_allele[2] = count_allele[2] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "GG") {
								count_geno[8] = count_geno[8] + 1;
								count_allele[3] = count_allele[3] + 2;
							} else if (best_geno[ig] == "GT") {
								count_geno[9] = count_geno[9] + 1;
								count_allele[3] = count_allele[3] + 1;
								count_allele[4] = count_allele[4] + 1;
							} else if (best_geno[ig] == "TT") {
								count_geno[10] = count_geno[10] + 1;
								count_allele[4] = count_allele[4] + 2;
							}
							null_error = (ind_coverage[ig] - read_count[ig][n1])/(double)ind_coverage[ig];
							if (n1 == 1) {
								null_maxll = read_count[ig][1]*log(1.0-null_error) + (ind_coverage[ig] - read_count[ig][1])*log(null_error/3.0);
							} else if (n1 == 2) {
								null_maxll = read_count[ig][2]*log(1.0-null_error) + (ind_coverage[ig] - read_count[ig][2])*log(null_error/3.0);
							} else if (n1 == 3) {
								null_maxll = read_count[ig][3]*log(1.0-null_error) + (ind_coverage[ig] - read_count[ig][3])*log(null_error/3.0);
							} else if (n1 == 4) {
								null_maxll = read_count[ig][4]*log(1.0-null_error) + (ind_coverage[ig] - read_count[ig][4])*log(null_error/3.0);
							}
							ind_llstat = 2.0*(maxll-null_maxll);
							if (best_geno[ig] == "AC") {
								if (ind_llstat > cv) {
									if (read_count[ig][1] >= read_count[ig][2]) {
										if (maxll - llg[1] <= cv) {
											test_geno = "AA";
										} else {
											test_geno = "AC";
										}
									} else {
										if (maxll - llg[5] <= cv) {
											test_geno = "CC";
										} else {
											test_geno = "AC";
										}
									}
								} else {
									test_geno = "AC";
								}
							} else if (best_geno[ig] == "AG") {
								if (ind_llstat > cv) {
                                                                        if (read_count[ig][1] >= read_count[ig][3]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "AA";
                                                                                } else {
                                                                                        test_geno = "AG";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[8] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "AG";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "AG";
                                                                }
							} else if (best_geno[ig] == "AT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][1] >= read_count[ig][4]) {
                                                                                if (maxll - llg[1] <= cv) {
                                                                                        test_geno = "AA";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[10] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "AT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "AT";
                                                                }
							} else if (best_geno[ig] == "CG") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][2] >= read_count[ig][3]) {
                                                                                if (maxll - llg[5] <= cv) {
                                                                                        test_geno = "CC";
                                                                                } else {
                                                                                        test_geno = "CG";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[8] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "CG";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "CG";
                                                                }
							} else if (best_geno[ig] == "CT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][2] >= read_count[ig][4]) {
                                                                                if (maxll - llg[5] <= cv) {
                                                                                        test_geno = "CC";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[10] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "CT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "CT";
                                                                }
							} else if (best_geno[ig] == "GT") {
                                                                if (ind_llstat > cv) {
                                                                        if (read_count[ig][3] >= read_count[ig][4]) {
                                                                                if (maxll - llg[8] <= cv) {
                                                                                        test_geno = "GG";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        } else {
                                                                                if (maxll - llg[10] <= cv) {
                                                                                        test_geno = "TT";
                                                                                } else {
                                                                                        test_geno = "GT";
                                                                                }
                                                                        }
                                                                } else {
                                                                        test_geno = "GT";
                                                                }
							} else {
								test_geno = best_geno[ig];
							}
							for (ng=0; ng<=1; ng++) {
								if (test_geno.substr(ng,1) != major_allele) {
									if (test_geno.substr(ng,1) == "A") {
										if (ind_llstat > max_ind_llstat[1]) {
											max_ind_llstat[1] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "C") {
										if (ind_llstat > max_ind_llstat[2]) {
											max_ind_llstat[2] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "G") {
										if (ind_llstat > max_ind_llstat[3]) {
											max_ind_llstat[3] = ind_llstat;
										}
									} else if (test_geno.substr(ng,1) == "T") {
										if (ind_llstat > max_ind_llstat[4]) {
											max_ind_llstat[4] = ind_llstat;
										}
									}
								}
							}
						} else {
							best_geno[ig] = "NA";
						}
					}
				}
				if (best_geno[ig] != "NA") {
					sum_best_error = sum_best_error + best_error;
				}
				/*
				if (max_ind_llstat[3] > cv) {
					printf("%d\n", ig);
				}
				*/
				/*
				if (ig == 1) {
					printf("%d\t%d\t%d\t%d\t%d\n", ig, read_count[ig][1], read_count[ig][2], read_count[ig][3], read_count[ig][4]);
					printf("best_error: %f\tllg[1]: %f\tllg[2]: %f\tllg[3]: %f\tbest_geno: %s\ttest:geno: %s\n", best_error, llg[1], llg[2], llg[3], best_geno[ig].c_str(), test_geno.c_str());
				}
				*/
  			}	// end of the loop over the individuals
			if (tot_count_geno > 0) {
				mean_best_error = sum_best_error/(double)tot_count_geno;
				num_alleles = 1;
				for(ag=1; ag<=4; ag++) {
					if (max_ind_llstat[ag] > cv) {
						num_alleles = num_alleles + 1;
					}
				}
				mean_best_error = sum_best_error/(double)tot_count_geno;
				// printf("num_alleles: %d\n", num_alleles);
  				// printf("scaffold: %d\tposition: %d\tsample size: %d\tallele 1 homo: %f\thetero: %f\tallele 2 homo: %f\n", scaffold, site, tot_count_geno, (double)count_geno[1]/tot_count_geno, (double)count_geno[2]/tot_count_geno, (double)count_geno[3]/tot_count_geno);
				if (count_allele[1] >= count_allele[2]) {
					a1 = 1;
					a2 = 2;
				} else {
					a1 = 2;
					a2 = 1;
				}
				if (count_allele[3] > count_allele[a1]) {
					a3 = a2;
					a2 = a1;
					a1 = 3;
				} else if (count_allele[3] > count_allele[a2]) {
					a3 = a2;
					a2 = 3;
				} else {
					a3 = 3;
				}
				if (count_allele[4] > count_allele[a1]) {
					a4 = a3;
					a3 = a2;
					a2 = a1;
					a1 = 4;
				} else if (count_allele[4] > count_allele[a2]) {
					a4 = a3;
					a3 = a2;
					a2 = 4;
				} else if (count_allele[4] > count_allele[a3]) {
					a4 = a3;
					a3 = 4;
				} else {
					a4 = 4;
				}
				if (num_alleles == 1) {
					// monomorphic
					fprintf(outstream, "%s\t%d\t%s\t%d\t%s\tNA\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1], count_allele[n2], mean_best_error);
					// printf("%s\t%d\t%s\t%d\t%s\tNA\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1], count_allele[n2], mean_best_error);
					for(ig=1;ig<nsample;ig++) {
						fprintf(outstream, "%s\t", best_geno[ig].c_str());
						// printf("%s\t", best_geno[ig].c_str());
					}
					fprintf(outstream, "%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
					// printf("%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
  				} else if (num_alleles == 2 ) {
					// Two alleles segregating
					if (count_allele[n1] >= count_allele[n2]) {
						fprintf(outstream, "%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), minor_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1], count_allele[n2], mean_best_error);
						// printf("%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), minor_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1], count_allele[n2], mean_best_error);
					} else {
						fprintf(outstream, "%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, minor_allele.c_str(), major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n2], count_allele[n1], mean_best_error);
						// printf("%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, minor_allele.c_str(), major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n2], count_allele[n1], mean_best_error);
					}
					for(ig=1;ig<nsample;ig++) {
						fprintf(outstream, "%s\t", best_geno[ig].c_str());
						// printf("%s\t", best_geno[ig].c_str());
					}
					fprintf(outstream, "%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
					// printf("%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
  				} else if (num_alleles > 2) {
					// More than two alleles segregating
					fprintf(outstream, "%s\t%d\t%s\t%d\tNA\tNA\t%d\t%d\tNA\tNA\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, pop_coverage, tot_count_geno, mean_best_error);
					// printf("%s\t%d\t%s\t%d\tNA\tNA\t%d\t%d\tNA\tNA\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, pop_coverage, tot_count_geno, mean_best_error);
					for(ig=1;ig<nsample;ig++) {
                                		fprintf(outstream, "%s\t", best_geno[ig].c_str());
                                        	// printf("%s\t", best_geno[ig].c_str());
                                	}
                                	fprintf(outstream, "%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
                                	// printf("%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n", best_geno[nsample].c_str(), max_ind_llstat[1], max_ind_llstat[2], max_ind_llstat[3], max_ind_llstat[4], count_allele[a1], count_allele[a2], count_allele[a3], count_allele[a4]);
				}
			} else {
				fprintf(outstream, "%s\t%d\t%s\tNA\tNA\tNA\t%d\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
				// printf("%s\t%d\t%s\tNA\tNA\tNA\t%d\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
				for(ig=1;ig<nsample;ig++) {
					fprintf(outstream, "%s\t", best_geno[ig].c_str());
					// printf("%s\t", best_geno[ig].c_str());
				}
				fprintf(outstream, "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best_geno[nsample].c_str());
				// printf("%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best_geno[nsample].c_str());
			} 
		} else if (pop_read[n1] == pop_coverage && pop_coverage > 0) {
			tot_count_geno = 0;
			count_allele[n1] = 0;
			num_alleles = 1;
			for(ig=1;ig<=nsample;ig++) {
				if (ind_coverage[ig] >= min_cov && ind_coverage[ig] <= max_cov) {
					tot_count_geno = tot_count_geno + 1;
					count_allele[n1] = count_allele[n1] + 2;	
					if (n1 == 1) {
                                		best_geno[ig] = "AA";
                        		} else if (n1 == 2) {
                                		best_geno[ig] = "CC";
                        		} else if (n1 == 3) {
                                		best_geno[ig] = "GG";
                        		} else if (n1 == 4) {
                               			best_geno[ig] = "TT";
                        		}
				}
			}
			fprintf(outstream, "%s\t%d\t%s\t%d\t%s\tNA\t%d\t%d\t%d\t0\t0\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1]);
			// printf("%s\t%d\t%s\t%d\t%s\tNA\t%d\t%d\t%d\t0\t0\t", scaffold.c_str(), site, ref_nuc.c_str(), num_alleles, major_allele.c_str(), pop_coverage, tot_count_geno, count_allele[n1]);
			for(ig=1;ig<nsample;ig++) {
				fprintf(outstream, "%s\t", best_geno[ig].c_str());
				// printf("%s\t", best_geno[ig].c_str());
			}
			fprintf(outstream, "%s\tNA\tNA\tNA\tNA\t%d\t0\t0\t0\n", best_geno[nsample].c_str(), count_allele[n1]);
			// printf("%s\tNA\tNA\tNA\tNA\t%d\t0\t0\t0\n", best_geno[nsample].c_str(), count_allele[n1]);
		} else if (pop_coverage == 0) {
			fprintf(outstream, "%s\t%d\t%s\tNA\tNA\tNA\t%d\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
			// printf("%s\t%d\t%s\tNA\tNA\tNA\t%d\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
			for(ig=1;ig<nsample;ig++) {
                                fprintf(outstream, "%s\t", best_geno[ig].c_str());
                                // printf("%s\t", best_geno[ig].c_str());
                        }
			fprintf(outstream, "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best_geno[nsample].c_str());
			// printf("%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best_geno[nsample].c_str());
		}				
	}
  	return 0;  
}
