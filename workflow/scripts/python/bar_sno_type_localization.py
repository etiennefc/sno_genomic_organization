#!/usr/bin/python3
import pandas as pd
import functions as ft
import matplotlib.pyplot as plt
import collections as coll

# Load dfs
sno_type_df = pd.read_csv(snakemake.input.sno_type, sep='\t')
sno_type_df = sno_type_df.drop_duplicates(subset='gene_id')
sno_types = sorted(list(pd.unique(sno_type_df['snoRNA_type'])), reverse=True)
dfs = []
for p in snakemake.input.sno_localization:
    if 'giardia_lamblia' not in p:
        df = pd.read_csv(p, sep='\t')
        species = p.split('sno_cluster_')[1].replace('.tsv', '')
        df['species'] = species
        dfs.append(df)
location_df = pd.concat(dfs)

dfs2 = []
for p in snakemake.input.sno_pos:
    if 'giardia_lamblia' not in p:
        df = pd.read_csv(p, sep='\t')
        species = p.split('sno_localization_')[1].replace('.tsv', '')
        df['species'] = species
        dfs2.append(df)
sno_pos_df = pd.concat(dfs2)


# Define variables
animals = ["homo_sapiens", "macaca_mulatta", "mus_musculus", "rattus_norvegicus", "ornithorhynchus_anatinus", 
              "gallus_gallus", "danio_rerio", "xenopus_tropicalis", "caenorhabditis_elegans", "drosophila_melanogaster"]
plants = ["arabidopsis_thaliana", "oryza_sativa", "triticum_aestivum"]
fungi = ["schizosaccharomyces_pombe", "neurospora_crassa", "candida_albicans", "saccharomyces_cerevisiae"]
protists = ["dictyostelium_discoideum", "giardia_lamblia", "tetrahymena_thermophila"]
location_type = list(pd.unique(location_df['localization']))
colors_location_type_dict = snakemake.params.colors_location_type
colors_location_type = colors_location_type_dict.values()

# Drop MRP snoRNA in S. cerevisiae
location_df = location_df[location_df['gene_id'] != 'NME1']
sno_pos_df = sno_pos_df[sno_pos_df['gene_id'] != 'NME1']

# Get the number of snoRNA per specie and order it in decreasing order per eukaryotic kingdom
animal_df = location_df[location_df['species'].isin(animals)]
plant_df = location_df[location_df['species'].isin(plants)]
fungi_df = location_df[location_df['species'].isin(fungi)]
protist_df = location_df[location_df['species'].isin(protists)]
other_df = location_df[~location_df['species'].isin(animals)]

animal_sno_nb_dict = dict(coll.Counter(animal_df['species']).most_common())  # in decreasing order
other_sno_nb_dict = dict(coll.Counter(plant_df['species']).most_common())  # in decreasing order
for df in [fungi_df, protist_df]:
    temp_dict = dict(coll.Counter(df['species']).most_common())
    for k,v in temp_dict.items():
        if k not in other_sno_nb_dict.keys():
            other_sno_nb_dict[k] = v
other_sno_nb_dict['giardia_lamblia'] = 0  # add the fact that no snoRNAs are annotated in that species

# Count the number of snoRNA per location type per animal species
count_list_animal = []
for sp in animal_sno_nb_dict.keys():
    temp_l = []
    temp_df = animal_df[animal_df['species'] == sp]
    for loca in location_type:
        temp_l.append(len(temp_df[temp_df['localization'] == loca]))
    count_list_animal.append(temp_l)
percent_animal = ft.percent_count(count_list_animal)

# Count the number of snoRNA per sno_type per location type per animal species
count_list_type_animal = []
for sp in animal_sno_nb_dict.keys():
    temp_l = []
    temp_df = animal_df[animal_df['species'] == sp]
    for loca in location_type:
        temp_loca = temp_df[temp_df['localization'] == loca]
        for s_type in sno_types:
            temp_s = sno_type_df[(sno_type_df['gene_id'].isin(list(temp_loca['gene_id']))) & (sno_type_df['snoRNA_type'] == s_type)]
            temp_l.append(len(temp_s))
    count_list_type_animal.append(temp_l)
percent_type_animal = ft.percent_count(count_list_type_animal)

# Create grouped bar chart for snoRNAs in animal species
ft.grouped_stacked_bar2(percent_type_animal, percent_animal, list(animal_sno_nb_dict.keys()), sno_types * len(location_type), 
                location_type, 'Animals', 'Species', 'Proportion of snoRNAs (%)', 
                ['lightgrey', 'grey', '#262626'] * len(location_type), list(colors_location_type), 0, 110, 
                [f'({i})' for i in animal_sno_nb_dict.values()], snakemake.output.bar_animals, custom_legend=True, width=0.35)


# Count the number of snoRNA per location type per other species
count_list_other = []
for sp in other_sno_nb_dict.keys():
    temp_l = []
    temp_df = other_df[other_df['species'] == sp]
    for loca in location_type:
        temp_l.append(len(temp_df[temp_df['localization'] == loca]))
    count_list_other.append(temp_l)
percent_other = ft.percent_count(count_list_other)

# Count the number of snoRNA per sno_type per location type per other species
count_list_type_other = []
for sp in other_sno_nb_dict.keys():
    temp_l = []
    temp_df = other_df[other_df['species'] == sp]
    for loca in location_type:
        temp_loca = temp_df[temp_df['localization'] == loca]
        for s_type in sno_types:
            temp_s = sno_type_df[(sno_type_df['gene_id'].isin(list(temp_loca['gene_id']))) & (sno_type_df['snoRNA_type'] == s_type)]
            temp_l.append(len(temp_s))
    count_list_type_other.append(temp_l)
percent_type_other = ft.percent_count(count_list_type_other)

# Create grouped bar chart for snoRNAs in other species
ft.grouped_stacked_bar2(percent_type_other, percent_other, list(other_sno_nb_dict.keys()), sno_types * len(location_type), 
                location_type, 'Plants, Fungi, Protists', 'Species', 'Proportion of snoRNAs (%)', 
                ['lightgrey', 'grey', '#262626'] * len(location_type), list(colors_location_type), 0, 110, 
                [f'({i})' for i in other_sno_nb_dict.values()], snakemake.output.bar_others, custom_legend=True, width=0.35)


# Format final_df 
final_df = sno_pos_df.merge(location_df[['gene_id', 'localization']], on='gene_id', how='left')
final_df = final_df.merge(sno_type_df[['gene_id', 'snoRNA_type']], on='gene_id', how='left')
final_df.loc[final_df['species'].isin(animals), 'eukaryotic_kingdom'] = 'Animals'
final_df.loc[final_df['species'].isin(plants), 'eukaryotic_kingdom'] = 'Plants'
final_df.loc[final_df['species'].isin(fungi), 'eukaryotic_kingdom'] = 'Fungi'
final_df.loc[final_df['species'].isin(protists), 'eukaryotic_kingdom'] = 'Protists'
final_df = final_df.sort_values(by=['eukaryotic_kingdom', 'species', 'gene_id'])


final_df = final_df[['gene_id', 'gene_name', 'species', 'eukaryotic_kingdom', 'localization', 
                    'snoRNA_type', 'gene_id_HG', 'gene_name_HG', 'gene_biotype_HG']]
print(final_df)
final_df.to_csv(snakemake.output.final_df, sep='\t', index=False)




