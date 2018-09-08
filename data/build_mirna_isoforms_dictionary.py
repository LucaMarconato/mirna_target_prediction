#!/usr/bin/env python3.7
import os
import re
from collections import OrderedDict
from typing import Dict, Tuple

# prefixes = set()
# the interpreter gives me an error on the following commented line
# mirna_isoforms: OrderedDict[Tuple[str], str] = OrderedDict()
mirna_isoforms = OrderedDict()
with open('raw/mirna.dat', 'r') as infile:
	accession_id_just_found = False
	for line in infile:
		line = line.rstrip('\n')
		prefix = line[0:2]
		# prefixes.add(prefix)
		if prefix == 'ID':
			r = re.compile(r'^ID\s+([a-zA-Z0-9\-]*?)\s+.*$')
			if r.match(line) == None:
				print(f'error: regexp does not match "{line}"')
				os._exit(1)
			ambiguous_mirna_id = r.sub(r'\1', line)
		elif prefix == 'FT':
			r = re.compile(r'^FT\s+/accession="(.*?)"$')
			if r.match(line) != None:
				accession_id_just_found = True
				accession_id = r.sub(r'\1', line)
			elif accession_id_just_found:
				accession_id_just_found = False
				r = re.compile(r'^FT\s+/product="(.*?)"$')
				if r.match(line) == None:
					print(f'error: regexp does not match "{line}"')
					os._exit(1)
				disambiguated_mirna_id = r.sub(r'\1', line)
				if(ambiguous_mirna_id[0:4] == 'hsa-'):
					mirna_isoforms[(ambiguous_mirna_id, accession_id)] = disambiguated_mirna_id
# print(prefixes)
# we sort lexicographically by the first and the second element of the key
# so x[0][0] is the ambiguous_mirna_id and x[0][1] is the accession_id
mirna_isoforms = sorted(mirna_isoforms.items(), key = lambda x: (x[0][0], x[0][1]))
with open('processed/mirna_isoforms_dictionary.tsv', 'w') as outfile:
	s = 'ambiguous_mirna_id\taccession_id\tdisambiguated_mirna_id\n'
	for key, value in mirna_isoforms:
		s += f'{key[0]}\t{key[1]}\t{value}\n'
	outfile.write(s)
print('done')

	
