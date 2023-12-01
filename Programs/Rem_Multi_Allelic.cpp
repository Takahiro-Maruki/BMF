// Updated on 09/19/16

/* 
	Program Rem_Multi_Allelic.cpp to remove sites containing
	more than two alleles from the population-level pro file.

	Inputs:
	1. Population-level pro file
	2. List of sites containing more than two alleles
*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* pro_file_name = {"pro.txt"};
	const char* ma_file_name = {"List_MA.txt"};
	const char* out_file_name = {"Out_RMA.txt"};
	int print_help = 0;

	FILE *outstream;	// pointer to the output file
	
	int argz = 1;	// argument counter

	// Read the specified setting
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-pf") == 0) {
			pro_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-mf") == 0) {
			ma_file_name = argv[++argz];
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
		fprintf(stderr, "       -pf <s>: specify the name of the pro file\n");
		fprintf(stderr, "       -mf <s>: specify the name of the file of multiple alleles\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		exit(1);
	}
	
	// Open the output file
	outstream = fopen(out_file_name, "w");
	if (outstream == NULL ) {	// Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file_name);
		exit(1);
	}

	ifstream pro_inputFile(pro_file_name);		// Try to open the input pro file
	if ( !pro_inputFile.is_open() ) {		// Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", pro_file_name); 
		exit(1);
	}
	
	// Read the header of the pro file
	string pro_line; // String buffer
	string h_scaf, h_site, h_ref_nuc; // Stores header names 
	getline(pro_inputFile, pro_line);
	istringstream pro_ss(pro_line);
	pro_ss >> h_scaf >> h_site >> h_ref_nuc;
	string str; // Temporarily stores each individual ID
	vector <string> id_ind; // Stores individual IDs.
	id_ind.clear();
	while (true) {
		pro_ss >> str;
		id_ind.push_back(str);
		if ( pro_ss.eof() ) {
			break;
		}
	}
	int nsample = id_ind.size();
	printf("%d individuals in the pro file\n", nsample);
	fprintf(outstream, "%s\n", pro_line.c_str());
	// printf("%s\n", pro_line.c_str());

	ifstream ma_inputFile(ma_file_name);            // Try to open the input file of multiple alleles
        if ( !ma_inputFile.is_open() ) {                // Exit on failure
                fprintf(stderr, "Cannot open %s for reading.\n", ma_file_name);
                exit(1);
        }

	// Read the header of the input file of multiple alleles
        string ma_line; // String buffer
        getline(ma_inputFile, ma_line);
        // printf("%s\n", ma_line.c_str());

	// Remove sites with more than two alleles from the pro file
	string ma_scaf;
	int ma_site;
	int matched;
	string pro_scaf;
	int pro_site;
	string ref_nuc;
	int ig;		// individual counter

	while ( getline(ma_inputFile, ma_line) ) {
		istringstream ma_ss(ma_line);
		ma_ss >> ma_scaf >> ma_site;
		matched = 0;
		while (matched == 0) {
			getline(pro_inputFile, pro_line);
			istringstream pro_ss(pro_line);
			pro_ss >> pro_scaf >> pro_site >> ref_nuc;
			if (pro_scaf == ma_scaf && pro_site == ma_site) {
				fprintf(outstream, "%s\t%d\t%s\t", pro_scaf.c_str(), pro_site, ref_nuc.c_str());
				// printf("%s\t%d\t%s\t", pro_scaf.c_str(), pro_site, ref_nuc.c_str());
				for(ig=0; ig<nsample-1; ig++) {
					fprintf(outstream, "0/0/0/0\t");
					// printf("0/0/0/0\t");
				}
				fprintf(outstream, "0/0/0/0\n");
				// printf("0/0/0/0\n");
				matched = 1;
			} else {
				fprintf(outstream, "%s\n", pro_line.c_str());
				// printf("%s\n", pro_line.c_str());
			}
		}
	}
	while (getline(pro_inputFile, pro_line) ) {
		fprintf(outstream, "%s\n", pro_line.c_str());
		// printf("%s\n", pro_line.c_str());
	}
	
	// Close the input and output files
	pro_inputFile.close();
	ma_inputFile.close();
	fclose(outstream);
	
	return(0);
}
