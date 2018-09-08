#!/usr/bin/env python3.7
import csv
import json
from collections import OrderedDict # to dump the json file respecting the order

l = list()
with open('raw/UTR_Sequences.txt', 'r') as infile:
	content = csv.reader(infile, delimiter='\t', quoting=csv.QUOTE_NONE)
	header = True
	
	# print(f'there are {len(content)} rows')
	processed = -1
	human_found = 0
	for row in content:
		processed += 1
		if (processed % 100000) == 0:
			print(f'{processed}/2300000 rows processed')
		if header:
			header = False
			continue
		e = dict()
		species_id = row[3]
		if species_id != '9606':
			continue		
		e['refseq_id'] = row[0]
		e['gene_id'] = row[1]
		e['gene_symbol'] = row[2]
		if e['refseq_id'] == 'ENST00000219343.6':
			print(f'found ENST00000219343.6, {processed} processed, {human_found} human_found')
		human_found += 1
		e['utr'] = row[4].replace('-','')		
		l.append(e)
		if (human_found % 1000) == 0:
			print(f'found {human_found} human utrs')

j = dict()		
j['list'] = l
with open('processed/human_utrs.json','w') as outfile:
	s = json.dumps(j, indent=4, sort_keys=False)
	outfile.write(s)

# import pdb; pdb.set_trace()
