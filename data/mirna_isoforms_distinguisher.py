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
		len_before = len(df)
		# print(f'len(df) = {len(df)} (before disambiguating)')
		global df_excluded
		df_excluded = df.iloc[0:0]
		for index, row in df.iterrows():
			mirna_id = index[0]
			mirna_region = index[1]
			split = mirna_region.split(',')
			if len(split) != 2:
				df_excluded = df_excluded.append(row)
				df.drop(index, inplace = True)
			else:
				if split[0] != 'mature':
					print('found mirna with accession number but not mature')
					df_excluded = df_excluded.append(row)
				else:
					df.rename({row.name[1]: split[1]}, axis = 'index', inplace = True, level = 1)
		sum_len_after = len(df) + len(df_excluded)
		# print(f'len(df) = {len(df)}, len(df_excluded) = {len(df_excluded)} (after disambiguating)')
		if len_before != sum_len_after:
			print(f'error: len_before = {len_before}, sum_len_after = {sum_len_after}')
			os._exit(1)
		df_isoform_dictionary = pd.read_csv('processed/mirna_isoforms_dictionary.tsv', delimiter = '\t')		
		key_checker = df_isoform_dictionary.groupby(['ambiguous_mirna_id', 'accession_id']).count()
		key_ambiguities_count = len(key_checker.loc[key_checker['disambiguated_mirna_id'] > 1, :])
		if key_ambiguities_count > 0:
			print(f'error: key_ambiguities_count = {key_ambiguities_count}')
			os._exit(1)
		df.reset_index(inplace = True)
		df.rename(columns = {'miRNA_ID': 'ambiguous_mirna_id', 'miRNA_region': 'accession_id'}, inplace = True)
		df_merged = pd.merge(df, df_isoform_dictionary, how = 'left', on = ['ambiguous_mirna_id', 'accession_id'])
		global df1
		df1 = df_merged
		
	else:
		print(f'skipping {tissue} mirnas')

distinguish_mirna_isoforms('normal')
distinguish_mirna_isoforms('tumor')
print('done')
