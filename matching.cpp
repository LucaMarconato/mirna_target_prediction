#include "matching.hpp"

#include <iostream>
#include <cstdio>
#include <cstring>

void valid_length(char * str0, char * str1, int length)
{
    int len0 = strlen(str0);
    int len1 = strlen(str1);
    if(len0 < length || len1 < length) {
        std::cerr << "error: len0 = " << len0 << ", len1 = " << len1 << "\n";
        exit(1);
    }
}

int dna_rna_mismatch(char * dna, char * rna, int length, bool watson_crick)
{
    valid_length(dna, rna, length);
    for(int i = 0; i < length; i++) {
        if(!(
             (dna[i] == 'A' && rna[i] == 'U') || 
             (dna[i] == 'C' && rna[i] == 'G') ||
             (dna[i] == 'G' && rna[i] == 'C') ||
             (dna[i] == 'T' && rna[i] == 'A') ||
             (!watson_crick && dna[i] == 'G' && rna[i] == 'U') ||
             (!watson_crick && dna[i] == 'T' && rna[i] == 'G')
             )
           ) {
            return i;
        }
    }
    return -1;
}

int dna_dna_mismatch(char * dna0, char * dna1, int length)
{
    valid_length(dna0, dna1, length);
    for(int i = 0; i < length; i++) {
        if(!(
             (dna0[i] == 'A' && dna1[i] == 'T') || 
             (dna0[i] == 'C' && dna1[i] == 'G') ||
             (dna0[i] == 'G' && dna1[i] == 'C') ||
             (dna0[i] == 'T' && dna1[i] == 'A')
             )
           ) {
            return i;
        }
    }
    return -1;
}

int rna_rna_mismatch(char * rna0, char * rna1, int length, bool watson_crick)
{
    valid_length(rna0, rna1, length);
    for(int i = 0; i < length; i++) {
        if(!(
             (rna0[i] == 'A' && rna1[i] == 'U') || 
             (rna0[i] == 'C' && rna1[i] == 'G') ||
             (rna0[i] == 'G' && rna1[i] == 'C') ||
             (rna0[i] == 'U' && rna1[i] == 'A') ||
             (!watson_crick && rna0[i] == 'G' && rna1[i] == 'U') ||
             (!watson_crick && rna0[i] == 'U' && rna1[i] == 'G')
             )
           ) {
            return i;
        }
    }
    return -1;
}
