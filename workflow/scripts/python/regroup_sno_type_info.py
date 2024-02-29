#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import collections as coll

# Load dfs
cols = ['gene_id', 'gene_name', 'snoRNA_type']
sno_type_other = pd.read_csv(snakemake.input.sno_type_other, sep='\t')
sno_type_human = pd.read_csv(snakemake.input.sno_type_human, sep='\t')
sno_type_human = sno_type_human.rename(columns={'snoDB ID': 'gene_id', 'Symbol': 'gene_name', 
                                                'Box Type': 'snoRNA_type'})
sno_type_human['gene_name'] = sno_type_human['gene_name'].fillna(sno_type_human['gene_id'])
sno_type_human = sno_type_human[cols]

sno_type_cerevisiae = pd.read_csv(snakemake.input.sno_type_cerevisiae, sep='\t')
sno_type_cerevisiae['gene_name'] = sno_type_cerevisiae['gene_id']
sno_type_cerevisiae['snoRNA_type'] = sno_type_cerevisiae['sno_type']
sno_type_cerevisiae = sno_type_cerevisiae[cols]

# Filter dfs to keep only snoRNAs
sno_type_human = sno_type_human[sno_type_human['snoRNA_type'].isin(['H/ACA', 'C/D', 'AluACA', 'unknown'])]
sno_type_human['snoRNA_type'] = sno_type_human['snoRNA_type'].replace('AluACA', 'H/ACA').replace('unknown', 'Unknown')
sno_type_cerevisiae = sno_type_cerevisiae[sno_type_cerevisiae['snoRNA_type'] != 'MRP']

sno_type_species = snakemake.input.rnacentral_sno_type
dfs = [sno_type_human, sno_type_cerevisiae]
other_species = ['drosophila_melanogaster', 'caenorhabditis_elegans', 
                'schizosaccharomyces_pombe']  # species snoRNA type found in Dieci et al. 2009 Genomics
for path in sno_type_species:
    df = pd.read_csv(path, sep='\t')
    for sp in other_species:
        if sp in path:
            if sp == 'schizosaccharomyces_pombe':
                df = df.merge(sno_type_other[['gene_id', 'snoRNA_type2']], how='left', on='gene_id')
                df['snoRNA_type'] = df['snoRNA_type'].replace('Unknown', None)
                df['snoRNA_type'] = df['snoRNA_type'].fillna(df['snoRNA_type2'])
                df['snoRNA_type'] = df['snoRNA_type'].fillna('Unknown')
                df = df.drop(columns='snoRNA_type2')
            if sp == 'caenorhabditis_elegans':
                df = df.merge(sno_type_other[['gene_name', 'snoRNA_type2']], how='left', on='gene_name')
                df['snoRNA_type'] = df['snoRNA_type'].replace('Unknown', None)
                df['snoRNA_type'] = df['snoRNA_type'].fillna(df['snoRNA_type2'])
                df['snoRNA_type'] = df['snoRNA_type'].fillna('Unknown')
                df = df.drop(columns='snoRNA_type2')
            if sp == 'drosophila_melanogaster':
                # Based on Dieci et al, gene name containing 'Me' are C/D (methylation), whereas 'Psi' are H/ACA (Psi for pseudouridylation)
                # Based on their dataset, 5 other snoRNAs with a slight variation in their name are present
                other_cd = ['snoRNA:291', 'snoRNA:684', 'snoRNA:185', 'snoRNA:229', 'snoRNA:83E4-5']
                df.loc[df['gene_name'].str.contains('Me|CD'), 'snoRNA_type'] = 'C/D'
                df.loc[df['gene_name'].str.contains('Psi|aca'), 'snoRNA_type'] = 'H/ACA'
                df.loc[df['gene_name'].isin(other_cd), 'snoRNA_type'] = 'C/D'
            break
    dfs.append(df)

final_df = pd.concat(dfs)
final_df.to_csv(snakemake.output.sno_type_all, sep='\t', index=False)