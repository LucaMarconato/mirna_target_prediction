#!/usr/bin/env python3.6
import os
import csv
import json
import pickle
from collections import OrderedDict # to dump the json file respecting the order
from typing import Dict, List, Tuple

class Graph:
	def __init__(self):
		self.mirna_site_arcs: Dict[Tuple[Mirna, Site], Mirna_site_arc] = dict()
		self.mirna_site_arcs_left: Dict[Mirna, List[Site]] = dict()
		self.mirna_site_arcs_right: Dict[Site, List[Mirna]] = dict()
		self.list_of_scores_without_corresponding_interactions: List[Tuple[Mirna, Site]] = list()
		self.list_of_scores_with_discording_conservation_information: List[Tuple[Mirna, Site]] = list()
		# self.site_gene_arcs_left: Dict[Site, List[Gene]] = dict()
		# self.site_gene_arcs_right: Dict[Gene, List[Site]] = dict()

	def add_mirna_site_arc(self, mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type, conserved):
		key0 = Mirna(mirna_family)
		key1 = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
		key = tuple((key0, key1))
		value = Mirna_site_arc(conserved)
		self.mirna_site_arcs[key] = value
		
		if not key0 in self.mirna_site_arcs_left:
			self.mirna_site_arcs_left[key0] = list()
		self.mirna_site_arcs_left[key0].append(key1)
		
		if not key1 in self.mirna_site_arcs_right:
			self.mirna_site_arcs_right[key1] = list()
		self.mirna_site_arcs_right[key1].append(key0)

		# if not key1 in self.site_gene_arcs_left:
		# 	self.site_gene_arcs_left[key1] = list()
		# self.site_gene_arcs_left[key1].append(gene)

	def add_context_scores_for_mirna_site_arc(self, mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type, context_score, weighted_context_score, conserved):
		key0 = Mirna(mirna_family)
		key1 = Site(gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
		key = tuple((key0, key1))
		if not key in self.mirna_site_arcs:
			# print(f'error: unable to find {key} into self.mirna_site_arcs for mirna {mirna_family}')
			# os._exit(1)
			self.add_mirna_site_arc(mirna_family, gene_id, gene_symbol, transcript_id, utr_start, utr_end, seed_match_type)
			self.list_of_scores_without_corresponding_interactions.append(key)
		# else:
		all_right = self.mirna_site_arcs[key].add_context_scores(context_score, weighted_context_score, conserved)
		if not all_right:
			self.list_of_scores_without_corresponding_interactions.append(key)

class Mirna:
	def __init__(self, mirna_family):
		self.mirna_family = mirna_family

	def __hash__(self):
		return hash(self.mirna_family)

	def __eq__(self, other):
		return (self.mirna_family == other.mirna_family)

	def __ne__(self, other):
		return not(self == other)

# class Gene:
# 	def __init__(self, gene_id, gene_symbol, transcript_id):
# 		self.gene_id = gene_id
# 		self.gene_symbol = gene_symbol
# 		self.transcript_id = transcript_id
		
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

			mirna_family = 'hsa-' + row[0]
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
			if (found_interactions % 100000) == 0:
				print(f'found_interactions = {found_interactions}')

		expected_number_of_interactions = len(g.mirna_site_arcs) - original_count_of_mirna_site_arcs
		if found_interactions != expected_number_of_interactions:
			print(f'error: found_interactions = {found_interactions}, expected_number_of_interactions = {expected_number_of_interactions})')
			os._exit(1)
	print(f'{len(g.mirna_site_arcs)} interactions, {len(g.mirna_site_arcs_left)} involved mirnas, {len(g.mirna_site_arcs_right)} involved sites')

def update_g_with_scores(file_containing_scores, conserved: bool):
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
			if (found_scores % 100000) == 0:
				print(f'found_scores = {found_scores}')
			if (len(g.list_of_scores_without_corresponding_interactions) % 100000) == 0:
				print(f'len(g.list_of_scores_without_corresponding_interactions) = {len(g.list_of_scores_without_corresponding_interactions)}')
		print(f'found_scores = {found_scores}')
		print(f'len(g.list_of_scores_without_corresponding_interactions) = {len(g.list_of_scores_without_corresponding_interactions)}')
			
if os.path.isfile('g.pkl'):
	with open('g.pkl', 'rb') as f:
		print('loading g.pkl... ', end = '', flush = True)			
		g, = pickle.load(f)
		print('DONE')
else:
	g = Graph()
	add_matchings_to_g('raw/Conserved_Family_Info.txt', True)
	add_matchings_to_g('raw/Nonconserved_Family_Info.txt', False)
	update_g_with_scores('raw/Conserved_Site_Context_Scores.txt', True)
	update_g_with_scores('raw/Nonconserved_Site_Context_Scores.txt', False)

if not os.path.isfile('g.pkl'):
	with open('g.pkl', 'wb') as f:
		print('saving g.pkl... ', end = '', flush = True)			
		pickle.dump([g], f)
		print('DONE')

# convert_matchings_to_json('raw/Conserved_Family_Info.txt_SMALL', 'raw/Conserved_Sites_Context_Scores.txt')

# l = list()
# 		e['refseq_id'] = row[0]
# 		e['gene_id'] = row[1]
# 		e['gene_symbol'] = row[2]
# 		if e['refseq_id'] == 'ENST00000219343.6':
# 			print(f'found ENST00000219343.6, {processed} processed, {human_found} human_found')
# 		human_found += 1
# 		e['utr'] = row[4].replace('-','')		
# 		l.append(e)
# 		if (human_found % 1000) == 0:
# 			print(f'found {human_found} human utrs')

# j = dict()		
# j['list'] = l
# with open('processed/human_utrs.json','w') as outfile:
# 	s = json.dumps(j, indent=4, sort_keys=False)
# 	outfile.write(s)

# # import pdb; pdb.set_trace()
