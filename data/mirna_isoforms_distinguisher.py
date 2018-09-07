#!/usr/bin/env python3.7
import sys
import os
import json
import pandas as pd

sys.argv.append('TCGA-CJ-4642')

if len(sys.argv) != 2:
	print('usage: mirna_isoform_distibguisher.py <patient_id>')
	os._exit(1)

patient_id = sys.argv[1]
patient_folder = 'patients/' + patient_id + '/'
if not os.path.isdir(patient_folder):
	print(f'patient {patient_id} not found')
	os._exit(1)

json_file = patient_folder + 'info.json'
with open(json_file, 'r') as infile:
	data = json.load(infile)

def distinguish_mirna_isoforms(tissue: str):
	j = data[f'{tissue}_mirnas_original_isoform_quantification']
	isoforms_filename = patient_folder + j['uuid'] + '/' + j['file']
	if os.path.isfile(isoforms_filename):
		global df
		df = pd.read_csv(isoforms_filename, delimiter = '\t')
		df = df.groupby(['miRNA_ID','miRNA_region']).sum()
		global df_excluded
		df_excluded = df.iloc[0:0]
		for index, row in df.iterrows():
			mirna_id = index[0]
			mirna_region = index[1]
			split = mirna_region.split(',')
			if len(split) != 2:
				df_excluded = df_excluded.append(row)
			else:
				print('TODO')
		print(f'to_process = {to_process}')
			
	else:
		print(f'skipping {tissue} mirnas')

distinguish_mirna_isoforms('normal')
distinguish_mirna_isoforms('tumor')
print('done')


