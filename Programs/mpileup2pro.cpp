// Updated on 05/08/22

/*
Program mpileup2pro.cpp to prepare a pro file of nucleotide-read counts in the
GFE format from an mpileup file of multiple individuals.

Inputs:
1. Reference nucleotide file.
2. List of individual IDs.
3. mpileup file of multiple individuals.

Output: pro file of nucleotide-read counts in the GFE format.

*/

#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
using namespace std;

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* ref_file = {"RefNuc.txt"};
	const char* list_id = {""};
	const char* mp_file = {""};
	const char* out_file_name = {"Out.pro"};
	int print_help = 0;

	int argz = 1;	// argument counter

	// Read the specified setting
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-ref") == 0) {
			ref_file = argv[++argz];
		} else if (strcmp(argv[argz], "-id") == 0) {
			list_id = argv[++argz];
		} else if (strcmp(argv[argz], "-mp") == 0) {
			mp_file = argv[++argz];
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
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
		fprintf(stderr, "       -ref <s>: specify the name of the reference file\n");
		fprintf(stderr, "       -id <s>: specify the name of the list file of individual IDs\n");
		fprintf(stderr, "       -mp <s>: specify the name of the mpileup file\n");
		fprintf(stderr, "       -out <s>: specify the name of the output file\n");
		exit(1);
	}

	// point to the output file
	FILE *outstream;

	ifstream i_inputFile(list_id);		// Try to open the list of individual IDs
	if ( !i_inputFile.is_open() ) {		// Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", list_id);
		exit(1);
	}

	string i_line;
	vector <string> id_ind; // Stores individual IDs.
	int nsample;	// number of individuals

	// Clear the vector
	id_ind.clear();
	while ( getline(i_inputFile, i_line) ) {
		id_ind.push_back(i_line);
	}
	nsample = id_ind.size();
	printf("%d individuals analyzed\n", nsample);

	// Close the list of individual IDs
	i_inputFile.close();

	// Open the output file
	outstream = fopen(out_file_name, "w");
	if (outstream == NULL ) {	// Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file_name);
		exit(1);
	}

	// Print the header of the output file
	int ig; // individul counter

	fprintf(outstream, "scaffold\tsite\tref_nuc\t");
	// printf("scaffold\tsite\tref_nuc\t");
	for(ig=0; ig<nsample-1; ig++) {
		fprintf(outstream, "%s\t", id_ind.at(ig).c_str());
		// printf("%s\t", id_ind.at(ig).c_str());
	}
	fprintf(outstream, "%s\n", id_ind.at(ig).c_str());
        // printf("%s\n", id_ind.at(ig).c_str());

	ifstream r_inputFile(ref_file);         // Try to open the reference file
        if ( !r_inputFile.is_open() ) {         // Exit on failure
                fprintf(stderr, "Cannot open %s for reading.\n", ref_file);
                exit(1);
        }

	string r_line; // String buffer

	// Read the header
        getline(r_inputFile, r_line);

	ifstream m_inputFile(mp_file);          // Try to open the mpileup file
	if ( !m_inputFile.is_open() ) {         // Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", mp_file);
		exit(1);
	}

	string m_line;	// String buffer
	string m_scaf;
	string scaffold, ref_nuc;
	int m_site, site;
	string m_ref_nuc;
	int read;	// one if the reference file needs to be read, zero, otherwise
	int nuc_read[5];	// an array of nucleotide-read counts in each individual
	int rg, ng;	// counters
	int dp; // depth
	string nr;	// nucleotide reads
	string bq;	// base quality
	int length_nr;	// number of characters in nr
	int indel;	// one if the letter in nr is involved in indels, zero, otherwise
	int find_size;	// one if the size of the intron needs to be found, zero, otherwise
	string s_size_indel;
	int size_indel;
	int nag;
	int end_ng_indel;

	while ( getline(m_inputFile, m_line) ) {	// Count nucleotide reads at each site in the reference sequence
		istringstream m_ss(m_line);
		m_ss >> m_scaf >> m_site >> m_ref_nuc;
		read = 1;
		while (read == 1) {
			getline(r_inputFile, r_line);
			istringstream r_ss(r_line);
			r_ss >> scaffold >> site >> ref_nuc;
			if (scaffold == m_scaf && site == m_site) {	// Site with data found
				read = 0;
				// Print out the reference info
				fprintf(outstream, "%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
				// printf("%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
				for(ig=0; ig<nsample; ig++) {
					// Set the initial conditions
					for (rg=1; rg<5; rg++) {
						nuc_read[rg] = 0;
					}
					m_ss >> dp >> nr >> bq;
					if (dp >= 1) {	// At least one nuceotide-read
						length_nr = nr.length();
						indel = 0;
						find_size = 0;
						for (ng=0; ng<length_nr; ng++) {
							if (nr.substr(ng,1) == "." || nr.substr(ng,1) == ",") {	// Reference match
								if (m_ref_nuc == "A") {
									nuc_read[1] = nuc_read[1] + 1;
								} else if (m_ref_nuc == "C") {
                                                                        nuc_read[2] = nuc_read[2] + 1;
								} else if (m_ref_nuc == "G") {
                                                                        nuc_read[3] = nuc_read[3] + 1;
								} else if (m_ref_nuc == "T") {
                                                                        nuc_read[4] = nuc_read[4] + 1;
								}
							} else {
								if (nr.substr(ng,1) == "+" || nr.substr(ng,1) == "-") {
									indel = 1;
									find_size = 1;
									s_size_indel = "";
									nag = 1;
									while (find_size == 1) {
										if ( 0 <= (nr.at(ng+nag) - '0') && (nr.at(ng+nag) - '0') <= 9 ) {
											s_size_indel = s_size_indel + nr.substr(ng+nag,1);
											nag = nag + 1;
										} else {
											find_size = 0;
										}
									}
									size_indel = atoi( s_size_indel.c_str() );
									end_ng_indel = ng + (nag - 1) + size_indel;
									// printf("%d bp indel in individual %d at site %d on %s: end_ng_indel %d\n", size_indel, ig+1, site, scaffold.c_str(), end_ng_indel);
									// printf("site: %d\tstart_ng_indel: %d\tsize_indel: %d\tend_ng_indel: %d\n", site, start_ng_indel, size_indel, end_ng_indel);
								}
								if (indel == 1 && ng > end_ng_indel) {
									indel = 0;
								}
								if (indel == 0) {
									if (ng == 0 || nr.substr(ng-1, 1) != "^") {
										if (nr.substr(ng,1) == "A" || nr.substr(ng,1) == "a") {
											nuc_read[1] = nuc_read[1] + 1;
										} else if (nr.substr(ng,1) == "C" || nr.substr(ng,1) == "c") {
                                                                        		nuc_read[2] = nuc_read[2] + 1;
                                                        			} else if (nr.substr(ng,1) == "G" || nr.substr(ng,1) == "g") {
                                                                        		nuc_read[3] = nuc_read[3] + 1;
                                                        			} else if (nr.substr(ng,1) == "T" || nr.substr(ng,1) == "t") {
                                                                        		nuc_read[4] = nuc_read[4] + 1;
                                                        			}
									}
								}
							}
						}
					}
					// Print out the nucleotide read quartet for the individual
					if (ig == nsample - 1) {
                                		fprintf(outstream, "%d/%d/%d/%d\n", nuc_read[1], nuc_read[2], nuc_read[3], nuc_read[4]);
                                		// printf("%d/%d/%d/%d\n", nuc_read[1], nuc_read[2], nuc_read[3], nuc_read[4]);
                        		} else {
                                		fprintf(outstream, "%d/%d/%d/%d\t", nuc_read[1], nuc_read[2], nuc_read[3], nuc_read[4]);
                                		// printf("%d/%d/%d/%d\t", nuc_read[1], nuc_read[2], nuc_read[3], nuc_read[4]);
                       			}
				}
			} else {
				// Print out the reference info
                                fprintf(outstream, "%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
                                // printf("%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
				// Print out the nucleotide read quartets
				for(ig=0; ig<nsample; ig++) {
					if (ig == nsample - 1) {
						fprintf(outstream, "0/0/0/0\n");
						// printf("0/0/0/0\n");
					} else {
						fprintf(outstream, "0/0/0/0\t");
                                                // printf("0/0/0/0\t");
					}
				}
			}
		}
	}	// End of the loop over mpileup lines
	while ( getline(r_inputFile, r_line) ) {	// Print out the remaining sites in the reference
		istringstream r_ss(r_line);
                r_ss >> scaffold >> site >> ref_nuc;
		// Print out the reference info
                fprintf(outstream, "%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
                // printf("%s\t%d\t%s\t", scaffold.c_str(), site, ref_nuc.c_str());
                // Print out the nucleotide read quartets
                for(ig=0; ig<nsample; ig++) {
                	if (ig == nsample - 1) {
                        	fprintf(outstream, "0/0/0/0\n");
                                // printf("0/0/0/0\n");
                        } else {
                        	fprintf(outstream, "0/0/0/0\t");
                                // printf("0/0/0/0\t");
                        }
                }
	}
	// Close the input and output files
	m_inputFile.close();
	r_inputFile.close();
	fclose(outstream);

	return(0);
}
