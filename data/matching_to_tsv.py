#!/usr/bin/env python3.6
import os
import csv
import json
import pickle
import re
import subprocess
from collections import OrderedDict # to dump the json file respecting the order
from typing import Dict, List, Tuple, Set

from matching_to_json_extract_data import show_g_statistics

class Graph:
	def __init__(self):
		self.mirna_site_arcs: Dict[Tuple[Mirna, Site], Mirna_site_arc] = dict()
		self.mirna_site_arcs_left: Dict[Mirna, List[Site]] = dict()
		self.mirna_site_arcs_right: Dict[Site, List[Mirna]] = dict()
		self.list_of_scores_without_corresponding_interactions: List[Tuple[Mirna, Site]] = list()
		self.list_of_scores_with_corresponding_interactions: List[Tuple[Mirna, Site]] = list()
		self.list_of_scores_with_discording_conservation_information: List[Tuple[Mirna, Site]] = list()

	def add_mirna_site_arc(self, mirna_family: str, gene_id: str, gene_symbol: str, transcript_id: str, utr_start: str, utr_end: str,
						   seed_match_type: str,
						   conserved: bool):
		key0 = Mirna(mirna_family)
		key1 = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
		key = tuple((key0, key1))
		value = Mirna_site_arc(conserved)
		if value in self.mirna_site_arcs:
			print('warning!')
		self.mirna_site_arcs[key] = value
		
		if not key0 in self.mirna_site_arcs_left:
			self.mirna_site_arcs_left[key0] = list()
		self.mirna_site_arcs_left[key0].append(key1)
		
		if not key1 in self.mirna_site_arcs_right:
			self.mirna_site_arcs_right[key1] = list()
		self.mirna_site_arcs_right[key1].append(key0)

	def add_context_scores_for_mirna_site_arc(self, mirna_family: str, gene_id: str, gene_symbol: str, transcript_id: str, utr_start: str, utr_end: str, seed_match_type: str, context_score: str, weighted_context_score: str, conserved: bool):
		key0 = Mirna(mirna_family)
		key1 = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
		key = tuple((key0, key1))
		if not key in self.mirna_site_arcs:
			self.add_mirna_site_arc(mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type, conserved)
			self.list_of_scores_without_corresponding_interactions.append(key)
		else:
			if self.mirna_site_arcs[key].context_scores_defined:
				print('error: context score already defined for this interaction')
				os._exit(1)
			self.list_of_scores_with_corresponding_interactions.append(key)			
		consistent_conservation_information = self.mirna_site_arcs[key].add_context_scores(context_score, weighted_context_score, conserved)
		if not consistent_conservation_information:
			self.list_of_scores_with_discording_conservation_information.append(key)

class Mirna:
	def __init__(self, mirna_family):
		self.mirna_family = mirna_family

	def __hash__(self):
		return hash(self.mirna_family)

	def __eq__(self, other):
		return (self.mirna_family == other.mirna_family)

	def __ne__(self, other):
		return not(self == other)

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

class Mirna_site_arc: 
	def __init__(self, conserved):
		self.context_scores_defined = False
		self.conserved = conserved

	def add_context_scores(self, context_score, weighted_context_score, conserved):
		self.context_scores_defined = True
		self.context_score = context_score
		self.weighted_context_score = weighted_context_score
		if self.conserved != conserved:
			# print(f'error: self.conserved = {self.conserved}, conserved = {conserved}')
			return False
		return True

def add_matchings_to_g(file_containing_matchings, conserved: bool):
	found_interactions = 0
	original_count_of_mirna_site_arcs = len(g.mirna_site_arcs)
	with open(file_containing_matchings, 'r') as infile:
		content = csv.reader(infile, delimiter = '\t', quoting = csv.QUOTE_NONE)
		header = True
		for row in content:
			if header:
				header = False
				continue			
			species_id = row[4]
			if species_id != '9606':
				continue

			mirna_family_temp = row[0]
			if re.search(r'miR-', mirna_family_temp) != None and re.search(r'let-', mirna_family_temp) != None:
				print(f'WARNING! {mirna_family_temp} not considered in the dataset, insert it manually')
				continue
			prefix = ''
			if re.match('miR-', mirna_family_temp) != None:
				mirna_family_temp = re.sub(r'^miR-', '', mirna_family_temp)
				prefix = 'miR-'
			if re.match('let-', mirna_family_temp) != None:
				mirna_family_temp = re.sub(r'^let-', '', mirna_family_temp)
				prefix = 'let-'

			if prefix == '':
				print(f'error: row[0] = {row[0]} does not start with miR- or let-')
				os._exit(1)
				
			real_mirnas = mirna_family_temp.split('/')

			for real_mirna in real_mirnas:
				mirna_family = 'hsa-' + prefix + real_mirna
				gene_id = row[1]
				gene_symbol = row[2]
				transcript_id = row[3]
				utr_start = row[5]
				utr_end = row[6]
				seed_match_type = row[9]
				if not seed_match_type in ["7mer-a1", "7mer-m8", "8mer"]:
					# 6 mer is not in the list, possibly beacuse of too many possible matchings
					# if we want to consider 6mer matchings we need to analyse the string matchings manually (a good idea is to use a suffix tree)
					continue
				pct = row[10] # unused for the moment

				g.add_mirna_site_arc(mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type, conserved)

				found_interactions += 1
				# if (found_interactions % 100000) == 0:
				# 	print(f'found_interactions = {found_interactions}')

		expected_number_of_interactions = len(g.mirna_site_arcs) - original_count_of_mirna_site_arcs
		if found_interactions != expected_number_of_interactions:
			print(f'error: found_interactions = {found_interactions}, expected_number_of_interactions = {expected_number_of_interactions})')
			os._exit(1)
	print(f'{len(g.mirna_site_arcs)} interactions, {len(g.mirna_site_arcs_left)} involved mirnas, {len(g.mirna_site_arcs_right)} involved sites')

def update_g_with_scores(file_containing_scores, conserved: bool):
	original_count_of_scores_without_corresponding_interaction = len(g.list_of_scores_without_corresponding_interactions)
	with open(file_containing_scores, 'r') as infile:
		found_scores = 0
		content = csv.reader(infile, delimiter = '\t', quoting = csv.QUOTE_NONE)
		header = True
		for row in content:
			if header:
				header = False
				continue
			species_id = row[3]
			if species_id != '9606':
				continue

			gene_id = row[0]
			gene_symbol = row[1]
			transcript_id = row[2]
			mirna_family = row[4]
			site_type_type = row[5]
			if site_type_type == '1':
				seed_match_type = "7mer-a1"
			elif site_type_type == '2':
				seed_match_type = "7mer-m8"
			elif site_type_type == '3':
				seed_match_type = "8mer"
			else:
				continue
			utr_start = row[6]
			utr_end = row[7]
			context_score = row[8]
			weighted_context_score = row[10]
			g.add_context_scores_for_mirna_site_arc(mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type, context_score, weighted_context_score, conserved)
			found_scores += 1
			# if (found_scores % 100000) == 0:
			# 	print(f'found_scores = {found_scores}')
			a = len(g.list_of_scores_without_corresponding_interactions) - original_count_of_scores_without_corresponding_interaction
			# if (a % 100000) == 0:
			# 	print(f'{a} scores has not a corresponding interaction')
		print(f'found_scores = {found_scores}')
		print(f'{a} scores has not a corresponding interaction')

def print_count_lines(filename: str):
	try:
		s = subprocess.check_output(f'cat {filename} | wc -l', shell = True)
		s = s.decode('utf-8').rstrip('\n')
		print(f'{filename} has{s} lines')
	except subprocess.CalledProcessError:
		print('error')

def show_g_statistics(g):
	print('*showing statistics*')
	# remember that not all the lines refers to human miRNAs
	print_count_lines('raw/Conserved_Family_Info.txt')
	print_count_lines('raw/Conserved_Site_Context_Scores.txt')
	print(f'len(g.list_of_scores_without_corresponding_interactions) = {len(g.list_of_scores_without_corresponding_interactions)}')
	print(f'len(g.list_of_scores_with_corresponding_interactions) = {len(g.list_of_scores_with_corresponding_interactions)}')
	the_sum = len(g.list_of_scores_without_corresponding_interactions) + len(g.list_of_scores_with_corresponding_interactions)
	print(f'the sum is {the_sum}')
	print(f'len(g.mirna_site_arcs) = {len(g.mirna_site_arcs)}, len(g.mirna_site_arcs_left) = {len(g.mirna_site_arcs_left)}, len(g.mirna_site_arcs_right) = {len(g.mirna_site_arcs_right)}')
	j = 0
	i = 0
	for key in g.mirna_site_arcs.keys():
		i += 1
		value = g.mirna_site_arcs[key]
		if value.context_scores_defined:
			j += 1
	print(f'{j}/{i} rows has the score ({round(j/i,2)} %)')
		
def export_g_data(g):
	considered_mirnas: Set[Mirna] = set()
	considered_sites: Set[Site] = set()
	
	with open('processed/scored_interactions.tsv', 'w') as outfile:
		i = 0
		to_write = 'mirna_family\tgene_id\tgene_symbol\ttranscript_id\tutr_start\tutr_end\tseed_match_type\tcontext_score\tweighted_context_score\tconserved\n'
		for key in g.mirna_site_arcs.keys():
			value = g.mirna_site_arcs[key]
			if value.context_scores_defined:
				e0 = key[0]
				e1 = key[1]
				considered_mirnas.add(e0)
				considered_sites.add(e1)
				to_write += f'{e0.mirna_family}\t{e1.gene_id}\t{e1.gene_symbol}\t{e1.transcript_id}\t{e1.utr_start}\t{e1.utr_end}\t{e1.seed_match_type}\t{value.context_score}\t{value.weighted_context_score}\t{value.conserved}\n'
				i += 1
		outfile.write(to_write)
		print(f'{i} rows written')

	# scored_interactions_left: Dict[Mirna, List[Site]] = dict()
	# scored_interactions_right: Dict[Site, List[Mirna]] = dict()
	# with open('processed/scored_interactions_left.json', 'w') as outfile:
	# 	for key in g.mirna_site_arcs.keys():
	# 		value = g.mirna_site_arcs[key]
	# 		if value.context_scores_defined:
	# 			e0 = key[0]
	# 			e1 = key[1]
				
	with open('processed/mirnas_with_scored_interactions.tsv', 'w') as outfile:
		i = 0
		to_write = 'mirna_family\n'
		for mirna in considered_mirnas:
			to_write += f'{mirna.mirna_family}\n'
			i += 1
		outfile.write(to_write)
		print(f'{i} rows written')

	with open('processed/sites_with_scored_interactions.tsv', 'w') as outfile:
		i = 0
		to_write = 'gene_id\tgene_symbol\ttranscript_id\tutr_start\tutr_end\tseed_match_type\n'
		for site in considered_sites:
			to_write += f'{e1.gene_id}\t{e1.gene_symbol}\t{e1.transcript_id}\t{e1.utr_start}\t{e1.utr_end}\t{e1.seed_match_type}\n'
			i += 1
		outfile.write(to_write)
		print(f'{i} rows written')
			
if os.path.isfile('g.pkl'):
	with open('g.pkl', 'rb') as f:
		print('loading g.pkl... ', end = '', flush = True)
		g, = pickle.load(f)
		print('DONE')
else:
	g = Graph()
	
	# add_matchings_to_g('raw/Conserved_Family_Info.txt_SMALL', True)
	# update_g_with_scores('raw/Conserved_Site_Context_Scores.txt_SMALL', True)
	
	add_matchings_to_g('raw/Conserved_Family_Info.txt', True)
	add_matchings_to_g('raw/Nonconserved_Family_Info.txt', False)
	update_g_with_scores('raw/Conserved_Site_Context_Scores.txt', True)
	update_g_with_scores('raw/Nonconserved_Site_Context_Scores.txt', False)

show_g_statistics(g)
export_g_data(g)
	
# if not os.path.isfile('g.pkl'):
# 	with open('g.pkl', 'wb') as f:
# 		print('saving g.pkl... ', end = '', flush = True)			
# 		pickle.dump([g], f)
# 		print('DONE')
