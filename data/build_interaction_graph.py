#!/usr/bin/env python3.6
import os
import csv
from typing import Dict, List, Tuple

from mirna import Mirna
from gene import Gene, Site

mirna_id_dictionary_left: Dict[Mirna, int] = dict()
mirna_id_dictionary_right: Dict[int, Mirna] = dict()
site_id_dictionary_left: Dict[Site, int] = dict()
site_id_dictionary_right: Dict[int, Site] = dict()
gene_id_dictionary_left: Dict[Gene, int] = dict()
gene_id_dictionary_right: Dict[int, Gene] = dict()

gene_id_site_id_arcs: Dict[int, List[int]] = dict()
site_id_gene_id_arcs: Dict[int, int] = dict()

class Mirna_site_arc:
	def __init__(self, mirna_id, site_id, context_score, weighted_context_score, conserved):
		self.mirna_id = mirna_id
		self.site_id = site_id
		self.context_score = context_score
		self.weighted_context_score = weighted_context_score
		self.conserved = conserved

class Matchings_graph:
	# mirna_site_arcs: (mirna_id, site_id) -> mirna_site_arc
	mirna_site_arcs: Dict[Tuple[int, int], Mirna_site_arc] = dict()
	# mirna_site_arcs_left: mirna_id -> List[site_id]
	mirna_site_arcs_left: Dict[int, List[int]] = dict()
	# mirna_site_arcs_right: site_id -> List[mirna_id]
	mirna_site_arcs_right: Dict[int, List[int]] = dict()
	# mirna_gene_arcs: (mirna_id, gene_id) -> List[site_id]
	mirna_gene_arcs: Dict[Tuple[int, int], List[int]]
	# mirna_gene_arcs_left: mirna_id -> List[gene_id]
	mirna_gene_arcs_left: Dict[int, List[int]] = dict()
	# mirna_gene_arcs_right: gene_id -> List[mirna_id]
	mirna_gene_arcs_right: Dict[int, List[int]] = dict()
	
	@classmethod
	def add_arc(cls, mirna_id, site_id, gene_id, mirna_site_arc):
		if (mirna_id, site_id) in cls.mirna_site_arcs:
			print(f'error: (mirna_id, site_id) = {(mirna_id, site_id)} is already in cls.mirna_site_arcs')
			os._exit(1)
		if mirna_id in cls.mirna_site_arcs_left and site_id in cls.mirna_site_arcs_left[mirna_id]:
			print(f'error: cls.mirna_site_arcs_left[{mirna_id}] already contains {site_id}')
			os._exit(1)
		if site_id in cls.mirna_site_arcs_right and mirna_id in cls.right_left[site_id]:
			print(f'error: cls.mirna_site_arcs_right[{site_id}] already contains {mirna_id}')
			os._exit(1)
		cls.mirna_site_arcs[(mirna_id, site_id)] = mirna_site_arc
		if not mirna_id in cls.mirna_site_arcs_left:
			cls.mirna_site_arcs_left[mirna_id] = list()
		cls.mirna_site_arcs_left[mirna_id].append(site_id)
		if not site_id in cls.mirna_site_arcs_right:
			cls.mirna_site_arcs_right[site_id] = list()
		cls.mirna_site_arcs_right[site_id].append(mirna_id)

		if not (mirna_id, gene_id) in cls.mirna_gene_arcs:
			cls.mirna_gene_arcs[(mirna_id, gene_id)] = list()
		if site_id in cls.mirna_gene_arcs[(mirna_id, gene_id)]:
			print(f'error: site_id = {site_id} already in cls.mirna_gene_arcs[({mirna_id}, {gene_id})]')
			os._exit(1)
		cls.mirna_gene_arcs[(mirna_id, gene_id)].append(site_id)

def build_mirna_id_dictionaries():
	with open('processed/mirnas_with_scored_interactions.tsv', 'r') as infile:
		content = csv.reader(infile, delimiter = '\t', quoting = csv.QUOTE_NONE)
		header = content.__next__()
		i = 0
		for row in content:
			mirna_family = row[0]
			if mirna_family in mirna_id_dictionary_left:
				print(f'error: {mirna_family} is already associated to {mirna_id_dictionary_left[mirna]} in mirna_id_dictionary_left')
				os._exit(1)
			mirna = Mirna(mirna_family)
			mirna_id_dictionary_left[mirna] = i
			mirna_id_dictionary_right[i] = mirna
			i += 1

def build_site_id_dictionaries():
	with open('processed/sites_with_scored_interactions.tsv', 'r') as infile:
		content = csv.reader(infile, delimiter = '\t', quoting = csv.QUOTE_NONE)
		header = content.__next__()
		i = 0
		j = 0
		for row in content:
			gene_id = row[0]
			gene_symbol = row[1]
			transcript_id = row[2]
			utr_start = row[3]
			utr_end = row[4]
			seed_match_type = row[5]
			site = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
			site_id_dictionary_left[site] = i
			site_id_dictionary_right[i] = site
			gene = Gene(gene_id, gene_symbol, transcript_id)
			i += 1
			if not gene in gene_id_dictionary_left:
				gene_id_dictionary_left[gene] = j
				gene_id_dictionary_right[j] = gene
				j += 1
				
			site_id = i - 1
			gene_id = gene_id_dictionary_left[gene]
			if not gene_id in gene_id_site_id_arcs:
				gene_id_site_id_arcs[gene_id] = list()
			gene_id_site_id_arcs[gene_id].append(site_id)
			if site_id in site_id_gene_id_arcs:
				print(f'error: site_id = {site_id} already in site_id_gene_id_arcs')
				os._exit(1)
			site_id_gene_id_arcs[site_id] = gene_id

def build_graph():
	with open('processed/scored_interactions.tsv', 'r') as infile:
		content = csv.reader(infile, delimiter = '\t', quoting = csv.QUOTE_NONE)
		header = content.__next__()
		i = 0
		for row in content:
			mirna_family = row[0]
			gene_id = row[1]
			gene_symbol = row[2]
			transcript_id = row[3]
			utr_start = row[4]
			utr_end = row[5]
			seed_match_type = row[6]
			context_score = row[7]
			weighted_context_score = row[8]
			conserved = row[9]
			mirna = Mirna(mirna_family)
			site = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
			gene = Gene(gene_id, gene_symbol, transcript_id)
			mirna_id = mirna_id_dictionary_left[mirna]
			site_id = site_id_dictionary_left[site]
			gene_id = gene_id_dictionary_left[gene]
			arc = Mirna_site_arc(context_score, weighted_context_score, conserved)
			Matchings_graph.add_arc(mirna_id, site_id, gene_id, arc)
			
if __name__ == '__main__':
	build_mirna_id_dictionaries()
	print('mirna dictionaries built')
	build_site_id_dictionaries()
	print('site dictionaries built')
	build_graph()
	print('graph built')
