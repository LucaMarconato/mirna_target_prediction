#!/usr/bin/env python3.7

# I used the name mirna_site for this file since site it is a standard library
class Site:
	def __init__(self, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type):
		self.gene_id = gene_id
		self.gene_symbol = gene_symbol
		self.transcript_id = transcript_id
		self.utr_start = utr_start
		self.utr_end = utr_end
		self.seed_match_type = seed_match_type

	def __hash__(self):
		return hash((self.gene_id, self.gene_symbol, self.transcript_id, self.utr_start, self.utr_end, self.seed_match_type))

	def __eq__(self, other):
		return ((self.gene_id, self.gene_symbol, self.transcript_id, self.utr_start, self.utr_end, self.seed_match_type) == (other.gene_id, other.gene_symbol, other.transcript_id, other.utr_start, other.utr_end, other.seed_match_type))

	def __ne__(self, other):
		return not(self == other)

class Gene:
	def __init__(self, gene_id, gene_symbol, transcript_id):
		self.gene_id = gene_id
		# we use the following variable for determining the identity of a gene
		# it is possible that genes with different ids share the same gene symbol
		# this is rare and is probably an error (as noted in gene_names_analyzer.r)
		# I have not checked if there are genes with the same names but different versions (TODO)
		self.gene_id_without_version = self.gene_id.split(".")[0]
		self.gene_symbol = gene_symbol
		self.transcript_id = transcript_id
		
	def __hash__(self):
		return hash(self.gene_id_without_version)
		# return hash((self.gene_id, self.gene_symbol, self.transcript_id))

	def __eq__(self, other):
		return (self.gene_id_without_version == other.gene_id_without_version)
		# return ((self.gene_id, self.gene_symbol, self.transcript_id) == (other.gene_id, other.gene_symbol, other.transcript_id))

	def __ne__(self, other):
		return not(self == other)
