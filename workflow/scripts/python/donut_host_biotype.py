#!/usr/bin/python3
import pandas as pd
import functions as ft
import matplotlib.pyplot as plt
import collections as coll

# Load dfs
dfs = []
for p in snakemake.input.sno_localization:
    if 'giardia_lamblia' not in p:
        df = pd.read_csv(p, sep='\t')
        species = p.split('sno_localization_')[1].replace('.tsv', '')
        df['species'] = species
        dfs.append(df)
location_df = pd.concat(dfs)

# Define variables
animals = ["homo_sapiens", "macaca_mulatta", "mus_musculus", "rattus_norvegicus", "ornithorhynchus_anatinus", 
              "gallus_gallus", "danio_rerio", "xenopus_tropicalis", "caenorhabditis_elegans", "drosophila_melanogaster"]
plants = ["arabidopsis_thaliana", "oryza_sativa", "triticum_aestivum"]
fungi = ["schizosaccharomyces_pombe", "neurospora_crassa", "candida_albicans", "saccharomyces_cerevisiae"]
protists = ["dictyostelium_discoideum", "giardia_lamblia", "tetrahymena_thermophila"]
hg_biotypes = ["protein_coding", "non_coding", "intergenic"]

# Drop MRP snoRNA in S. cerevisiae
location_df = location_df[location_df['gene_id'] != 'NME1']

# Get the number of snoRNA per specie and order it in decreasing order per eukaryotic kingdom
animal_df = location_df[location_df['species'].isin(animals)]
plant_df = location_df[location_df['species'].isin(plants)]
fungi_df = location_df[location_df['species'].isin(fungi)]
protist_df = location_df[location_df['species'].isin(protists)]
other_df = location_df[~location_df['species'].isin(animals)]  # everything that's not an animal

animal_sno_nb_dict = dict(coll.Counter(animal_df['species']).most_common())  # in decreasing order
other_sno_nb_dict = dict(coll.Counter(plant_df['species']).most_common())  # in decreasing order
for df in [fungi_df, protist_df]:
    temp_dict = dict(coll.Counter(df['species']).most_common())
    for k,v in temp_dict.items():
        if k not in other_sno_nb_dict.keys():
            other_sno_nb_dict[k] = v
other_sno_nb_dict['giardia_lamblia'] = 0  # add the fact that no snoRNAs are annotated in that species

print(location_df)

# Count the number of snoRNA per HG biotype per animal species and create donut charts
count_list_animal = []
for sp in animal_sno_nb_dict.keys():
    temp_l = []
    temp_df = animal_df[animal_df['species'] == sp]
    for biotype in hg_biotypes:
        temp_l.append(len(temp_df[temp_df['gene_biotype_HG'] == biotype]))
    count_list_animal.append(temp_l)
percent_animal = ft.percent_count(count_list_animal)

ft.pie_multiple(1, 10, percent_animal, hg_biotypes, ['#4d9221', '#b8e186', '#dfc27d'], 
            list(animal_sno_nb_dict.keys()), 'Animals', None, snakemake.output.donut_animals)

# Count the number of snoRNA per HG biotype per other species and create donut charts
count_list_other = []
for sp in other_sno_nb_dict.keys():
    temp_l = []
    temp_df = other_df[other_df['species'] == sp]
    for biotype in hg_biotypes:
        temp_l.append(len(temp_df[temp_df['gene_biotype_HG'] == biotype]))
    count_list_other.append(temp_l)
percent_other = ft.percent_count(count_list_other)
print(percent_other)
ft.pie_multiple(1, 10, percent_other, hg_biotypes, ['#4d9221', '#b8e186', '#dfc27d'], 
            list(other_sno_nb_dict.keys()), 'Plants, Fungi, Protists', None, snakemake.output.donut_others)






