import csv
import json
from collections import OrderedDict # to dump the json file respecting the order

l = list()
with open('raw/miR_Family_Info.txt', 'r') as infile:
	content = csv.reader(infile, delimiter='\t', quotechar='"')
	header = True
	for row in content:
		if header:
			header = False
			continue
		e = dict()
		e['id'] = row[0]
		e['s'] = row[4]
		l.append(e)

j = dict()
j['list'] = l
with open('processed/human_mirnas.json', 'w') as outfile:
	s = json.dumps(j, indent=4, sort_keys=False)
	outfile.write(s)
