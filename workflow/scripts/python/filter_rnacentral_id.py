#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import collections as coll

rnacentral_ids = snakemake.input.ids
species = str(snakemake.wildcards.species)
#taxon_id = snakemake.params.taxon_id
url_part = snakemake.params.url_part
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', 
                names=['chr', 'start', 'end', 'gene_id', 
                        'score', 'strand', 'gene_name'])

# For each sno, find if it has a RNAcentral id
if species == 'giardia_lamblia':  # no snoRNA annotated in that species
    sp.call(f'touch rna_central_search_{species}.tsv', shell=True)
else:
    sno_ids = '|'.join(list(sno_bed['gene_id']))
    sp.call(f'grep -E "{sno_ids}" {rnacentral_ids} > rna_central_search_{species}.tsv', shell=True)

# Load id conversion df and save it in right format
print(species)
print(sno_bed)
df = pd.read_csv(f'rna_central_search_{species}.tsv', sep='\t', 
                names=['rnacentral_id', 'source', 'gene_id', 'taxon_id', 'gene_biotype', 'gene_id2'])
df = df.drop_duplicates(subset='gene_id2')
if species in ['dictyostelium_discoideum']:
    # Species for which there is duplicate entries in Dyctibase and Ensembl_protists
    df = df.drop_duplicates(subset='rnacentral_id')
df[['gene_id', 'gene_id2', 'rnacentral_id', 'taxon_id', 'source']].to_csv(snakemake.output.id_table, sep='\t', index=False)
print(df)

# Create url table
if species == 'giardia_lamblia': # no snoRNA annotated in that species
    pd.DataFrame().to_csv(snakemake.output.urls, sep='\t', index=False)
else:
    with open(snakemake.output.urls, 'w') as f:
        for i, line in df.iterrows():
            f.write(f"{url_part}{line['rnacentral_id']}_{line['taxon_id']}\n")


sp.call(f'rm rna_central_search_{species}.tsv', shell=True)



