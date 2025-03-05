## Code used for filtering SortMeRNA blast-like outputs as described in the methods section from the article entitled 
## "Unveiling the rich and functionally diverse Southern Brazilian Highland Grasslands soil funga for promoting conservation".

## Contact: kelmermartinscunha@gmail.com


import pandas as pd
import numpy as np
import os
import re
import itertools


def load_blast_files(path_to_files):
    pwd = os.listdir(path_to_files)
    blast = {}
    for file in pwd:
        if '.blast' in file:
            blast[file] = pd.read_csv(path_to_files+'/'+file, sep = "\t",
                     names=['query', 'subject', 'percentageID', 'alignmentLength', 'mismatches', 'gapOpenings',
                           'qStart', 'qEnd', 'sStart', 'sEnd', 'eValue', 'bitScore','cigar', 'queryCoverage', 'strand'],
                     index_col=False)
    for i in blast:
            blast[i][['species', 'code', 'sh', 'uniteType', 'taxonomy']] = blast[i].subject.str.split("|", expand=True)
            blast[i][['kingdom', 'phylum', 'class', 'order',
            'family', 'genus', 'species']] = blast[i].taxonomy.str.split(";", expand=True)

    return blast


def filter_blast(blast, eValue, queryCoverage, length, id):
    for i in blast:
        blast[i] = (blast[i][(blast[i]['eValue'] <= eValue) & (blast[i]['queryCoverage'] >= queryCoverage) &
        (blast[i]['alignmentLength'] >= length) & (blast[i]['percentageID'] >= id)])
        blast[i][['species', 'code', 'sh', 'uniteType', 'taxonomy']] = blast[i].subject.str.split("|", expand=True)
        blast[i][['kingdom', 'phylum', 'class', 'order','family', 'genus',
        'species']] = blast[i].taxonomy.str.split(";", expand=True)
        blast[i] = blast[i][blast[i].phylum != 'p__unidentified'].reset_index()

    for key in blast:
        df = blast[key]
        dups = df[df.duplicated(subset='query', keep=False)]
        results  = []
        for query, group in dups.groupby('query'):
            ordered_group = group.sort_values(['percentageID', 'alignmentLength', 'bitScore', 'eValue', 'queryCoverage'],
            ascending=[False, False, False, True, False]).reset_index(drop=True)
            results.append(ordered_group.iloc[1:len(ordered_group)])
        merged = df.merge(pd.concat(results), on=list(df.columns), how='left', indicator=True)
        blast[key] = merged[merged['_merge'] == 'left_only']
        blast[key] = blast[key].drop(columns='_merge')


def output_results(blast):
    result_df = pd.DataFrame()
    for df_name, df in blast.items():
        counts_df = df['sh'].value_counts().reset_index()
        counts_df.columns = ['sh', df_name]
        if result_df.empty:
            result_df = counts_df
        else:
            result_df = pd.merge(result_df, counts_df, on='sh', how='outer')
    result_df = result_df.fillna(0)
    result_df['sh'] = result_df['sh'].str.replace(" ", "")
    result_df = result_df.set_index(result_df['sh'])
    result_df = result_df.drop('sh', axis=1)
    return result_df


def main():
    blast_result = load_blast_files("./blast_files")
    filter_blast(blast_result, 0.01, 90, 70, 90)

    result_csv = output_results(blast_result)
    ddOTU_table = result_csv.transpose()
    ddOTU_table.index.name = None
    ddOTU_table = ddOTU_table.rename_axis(None, axis=1)
    ddOTU_table.to_csv('ddOTU_table.csv')

if __name__ == "__main__":
    main()
