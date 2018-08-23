#ifndef MATCHING_H
#define MATCHING_H

/*
  args:
  dna
  	dna must be a string of dna
  rna
  	rna must be a string of rna
  length
	it is checked for a matching in the first length characters
  watson_crick
  	if true only watson_crick pairing are checked, otherwise also G-U and T-G pairings (TODO: does this make sense?)
  returned value:
  -1 if there is no mismatch, otherwise the index of the mistmatch
 */
int dna_rna_mismatch(char * dna, char * rna, int length, bool watson_crick);

/*
  similar to dna_rna_mismatch, only watson_crick pairings are considered
*/
int dna_dna_mismatch(char * dna0, char * dna1, int length);

/*
  similar to dna_rna_mismatch
 */
int rna_rna_mismatch(char * rna0, char * rna1, int length, bool watson_crick);

#endif // MATCHING_H
