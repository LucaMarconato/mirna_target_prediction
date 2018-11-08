import os
import json
import pandas as pd


def distinguish_mirna_isoforms(data_path: str):
    json_file = data_path + 'info.json'
    with open(json_file, 'r') as infile:
        data = json.load(infile)
    j = data[f'original_mirna_isoform_quantification']
    isoforms_filename = data_path + j['uuid'] + '/' + j['file']
    if not os.path.isfile(isoforms_filename):
        print(f'error: the path does not exists, data_path = {data_path}')
        os._exit(1)
    else:
        # global df
        df = pd.read_csv(isoforms_filename, delimiter='\t')
        df = df.groupby(['miRNA_ID', 'miRNA_region']).sum()
        len_before = len(df)
        total_reads = sum(df['read_count'])
        # print(f'len(df) = {len(df)} (before disambiguating)')
        # global df_excluded
        df_excluded = df.iloc[0:0]
        for index, row in df.iterrows():
            mirna_id = index[0]
            mirna_region = index[1]
            split = mirna_region.split(',')
            if len(split) != 2:
                df_excluded = df_excluded.append(row)
                df.drop(index, inplace=True)
            else:
                if split[0] != 'mature':
                    print('found mirna with accession number but not mature')
                    df_excluded = df_excluded.append(row)
                else:
                    df.rename({row.name[1]: split[1]}, axis='index', inplace=True, level=1)
        excluded_reads = sum(df_excluded['read_count'])
        print((f'{len(df_excluded)} mirna precursors excluded, '
               f'{int(excluded_reads)}/{total_reads} = {round(excluded_reads/total_reads,7)} excluded reads'))
        sum_len_after = len(df) + len(df_excluded)
        # print(f'len(df) = {len(df)}, len(df_excluded) = {len(df_excluded)} (after disambiguating)')
        if len_before != sum_len_after:
            print(f'error: len_before = {len_before}, sum_len_after = {sum_len_after}')
            os._exit(1)
        df_isoform_dictionary = pd.read_csv('../../processed/mirna_isoforms_dictionary.tsv', delimiter='\t')
        # check that there are no duplicate keys in the dictionary,
        # where the first two columns of the data frame are considered as a the key
        key_checker = df_isoform_dictionary.groupby(['ambiguous_mirna_id', 'accession_id']).count()
        key_ambiguities_count = len(key_checker.loc[key_checker['disambiguated_mirna_id'] > 1, :])
        if key_ambiguities_count > 0:
            print(f'error: key_ambiguities_count = {key_ambiguities_count}')
            os._exit(1)

        df.reset_index(inplace=True)
        df.rename(columns={'miRNA_ID': 'ambiguous_mirna_id', 'miRNA_region': 'accession_id'}, inplace=True)
        # global df_merged
        df_merged = pd.merge(df, df_isoform_dictionary, how='left', on=['ambiguous_mirna_id', 'accession_id'])
        # The first 4 columns of df_merged: ambiguous_mirna_id, accession_id, read_count, reads_per_million_miRNA_mapped
        # cannot contain nan values, the fifth (and last) column, disambiguated_mirna_id, can.
        # In such a case we print a warning and we set the disambiguated_mirna_id equal to the ambiguous_mirna_id

        # check that the first 4 columns do not contain nan values
        nan_count = len(df_merged) - df_merged.count()
        for i in range(0, 4):
            if nan_count[i] > 0:
                print(f'error: nan_count[i] = {nan_count[i]}')
                os._exit(1)
        df_not_disambiguated = df_merged[df_merged['disambiguated_mirna_id'].isnull()]
        if len(df_not_disambiguated) > 0:
            not_disambiguated_reads = sum(df_not_disambiguated['read_count'])
            print((f'unable to disambiguate {len(df_not_disambiguated)} mirna precursors, '
                   f'{not_disambiguated_reads}/{total_reads} = '
                   f'{round(not_disambiguated_reads/total_reads,7)} not disambiguated reads'))
            # print(f'for those rows we set disambiguated_mirna_id equal to ambiguous_mirna_id')
            # for index, row in df_not_disambiguated.iterrows():
            #     to_set = df_merged.at[index, 'ambiguous_mirna_id']
            #     df_merged.at[index, 'disambiguated_mirna_id'] = to_set
        # global df_to_export
        df_to_export = pd.DataFrame({'mirbase_id': df_merged['disambiguated_mirna_id'],
                                     'reads': df_merged['read_count']})
        df_to_export = df_to_export.groupby('mirbase_id').sum().reset_index()
        df_to_export['mirbase_id'] = df_to_export['mirbase_id'].str.lower()
        df_to_export['rpm'] = df_to_export['reads'] / df_to_export['reads'].sum() * 1000000
        filename = data_path.replace('raw/', 'processed/') + 'mirna_expression_profile.tsv'
        df_to_export.to_csv(filename, header=True, sep='\t', index=False)


data_paths = ['TCGA-CJ-4642/tumor/raw/']

for data_path in data_paths:
    if not os.path.isdir(data_path):
        print(f'data path {data_path} not found')
        os._exit(1)
    else:
        print(f'processing {data_path}')
        distinguish_mirna_isoforms(data_path)
print('done')
