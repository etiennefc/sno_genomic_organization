#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
from gtfparse import read_gtf
import collections as coll

species = str(snakemake.wildcards.species)
print(species)
intergenic_distance = snakemake.params.intergenic_distance
cols = ['seqname', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name']

if species in ["giardia_lamblia"]:  # no snoRNA annotated in these species
    pd.DataFrame().to_csv(snakemake.output.sno_cluster, sep='\t', index=False)
else:
    sno_pos = pd.read_csv(snakemake.input.sno_pos, sep='\t')
    sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', names=cols)
    HG_bed = pd.read_csv(snakemake.input.HG_bed, sep='\t', names=cols + ['gene_biotype'])

    # Separate intergenic and intronic snoRNAs
    intergenic_sno = sno_pos[sno_pos['gene_biotype_HG'] == 'intergenic']
    intronic_sno = sno_pos[sno_pos['gene_biotype_HG'] != 'intergenic']
    initial_intronic_sno = intronic_sno.copy()

    # Sort intergenic snoRNAs per position
    intergenic_sno = intergenic_sno.sort_values(by=['seqname', 'start']).reset_index(drop=True)

    ## Groupby snoRNAs per chr; define intergenic snoRNAs as in a cluster based on distance
    # Function to determine if two snoRNA genes are within 500nt
    def is_within(current_end, next_start, dist):
        return next_start - current_end <= dist

    cluster_number = 1
    for i, group in intergenic_sno.groupby('seqname'):  # groupby chr
        current_cluster_end = group.iloc[0]['end']
        cluster_gene_ids = []
        index = 0
        for j, row in group.iterrows():
            # test if sno end is at <=500nt from the next snoRNA start
            if is_within(current_cluster_end, row['start'], intergenic_distance):
                # Add gene id to cluster and update cluster end
                cluster_gene_ids.append(row['gene_id'])
                current_cluster_end = row['end']
                # Account for the last potential cluster of a chr
                if index == len(group) - 1:
                    if len(cluster_gene_ids) > 1:
                        intergenic_sno.loc[intergenic_sno['gene_id'].isin(cluster_gene_ids), 
                                        'cluster'] = f'intergenic_cluster_{cluster_number}'
            else:
                # Update cluster end to the next snoRNA to start a new potential cluster
                current_cluster_end = row['end']
                if len(cluster_gene_ids) > 1:
                    # Give to previous positive cluster a name
                    intergenic_sno.loc[intergenic_sno['gene_id'].isin(cluster_gene_ids), 
                                        'cluster'] = f'intergenic_cluster_{cluster_number}'
                    cluster_number += 1
                cluster_gene_ids = [row['gene_id']]
            index += 1
            

    # Define the rest of intergenic sno as not in a cluster
    if 'cluster' in intergenic_sno.columns:
        intergenic_sno['cluster'] = intergenic_sno['cluster'].fillna('intergenic')
        intergenic_sno.loc[intergenic_sno['cluster'].str.contains('intergenic_cluster_'), 'localization'] = 'intergenic_cluster'
        intergenic_sno.loc[intergenic_sno['cluster'] == 'intergenic', 'localization'] = 'mono_intergenic'
    else:
        intergenic_sno['cluster'] = 'intergenic'
        intergenic_sno['localization'] = 'mono_intergenic'




    # Create temporary gtf of only the host genes
    if len(intronic_sno) > 0:
        for i in pd.unique(intronic_sno.gene_id_HG):
            sp.call(f"""awk '$0 ~ "{i}"' {snakemake.input.gtf} >> hg_gtf_temp_{species}.gtf""", shell=True)

        # Select transcript -201 (or the 202, 203,...) as representative transcript of the host gene if it ecompasses all snoRNAs it embeds
        # Select its exons
        temp_hg_gtf = read_gtf(f'hg_gtf_temp_{species}.gtf')
        exon_hg = {}
        false_intronic_sno = [] # sno that overlap a host gene, but not an indivual transcript of that HG
        for k, group in temp_hg_gtf.groupby('gene_id'):
            group = group.reset_index(drop=True)
            # Add transcript_name column when missing
            if 'transcript_name' not in group.columns:
                group['transcript_name'] = group['transcript_id'] + '-201_added'
            if pd.unique(group['transcript_name'])[0] == '':
                group['transcript_name'] = group['transcript_id'] + '-201_added'

            trans_df = group[group['feature'] == 'transcript']
            exon_df = group[group['feature'] == 'exon']
            gene_id = pd.unique(exon_df['gene_id'])[0]
            
            transcripts = sorted([t for t in trans_df['transcript_name']])
            sno = intronic_sno[intronic_sno['gene_id_HG'] == gene_id]
            min_start, max_end = sno['start'].min(), sno['end'].max()  # start of most upstream sno and end of most downstream sno in that HG
            # Choose one transcript per HG based on name (201, 202, ...) and if it encompasses all of its sno
            temp_t = None
            for t in transcripts:
                if (list(trans_df[trans_df['transcript_name'] == t]['start'])[0] <= min_start) & (list(trans_df[trans_df['transcript_name'] == t]['end'])[0] >= max_end):
                    chosen_transcript = t
                    temp_t = chosen_transcript
                    exons = exon_df[exon_df['transcript_name'] == chosen_transcript]
                    exon_hg[gene_id] = exons
                    break
            # For long HG that have multiple snoRNAs encoded in different transcripts that
            # don't span all the gene length (ex: SNHG14 in human) (so a single transcript does not cover all snoRNAs)
            if temp_t == None:
                # Create a sorted dict of longest to shortest transcript per HG
                t_length_dict = {list(trans_df[trans_df['transcript_name'] == t]['transcript_name'])[0]: 
                                list(trans_df[trans_df['transcript_name'] == t]['end'])[0] - 
                                list(trans_df[trans_df['transcript_name'] == t]['start'])[0] for t in transcripts}
                sorted_t_length_dict = {i[0]:i[1] for i in sorted(t_length_dict.items(), key=lambda x:x[1], reverse=True)}
                linked_snoRNA = {}
                # Iterate for each sno over all trasncripts; take the first that contains an intron overlapping the sno;
                # if not, take the first transcript that encompasses entirely the sno, but an exon overlaps with the sno
                for snoRNA in sno['gene_id']:
                    potential_trans = []
                    for it, t2 in enumerate(sorted_t_length_dict.keys()):
                        if (snoRNA not in linked_snoRNA.keys()) & (len(linked_snoRNA.keys()) < len(sno)):
                            st = list(sno[sno['gene_id'] == snoRNA]['start'])[0]
                            nd = list(sno[sno['gene_id'] == snoRNA]['end'])[0]
                            if (list(trans_df[trans_df['transcript_name'] == t2]['start'])[0] <= st) \
                                    & (list(trans_df[trans_df['transcript_name'] == t2]['end'])[0] >= nd):
                                potential_trans.append(t2) # transcript encompasses the sno entirely
                                select_exons2 = exon_df[exon_df['transcript_name'] == t2]
                                exon_pairs2 = [[select_exons2.iloc[i], select_exons2.iloc[i+1]] for i in range(len(select_exons2) - 1)]
                                for pair in exon_pairs2:
                                    exon1, exon2 = pair[0], pair[1]
                                    # If an intron is found overlapping that sno, use that transcript as host gene
                                    if (exon1['end'] < st) & (exon2['start'] > nd):
                                        linked_snoRNA[snoRNA] = t2
                                        exons = exon_df[exon_df['transcript_name'] == t2]
                                        intronic_sno.loc[intronic_sno['gene_id'] == snoRNA, 'gene_id_HG'] = t2 + '_mod'
                                        if t2 not in exon_hg.keys():
                                            exon_hg[t2 + '_mod'] = exons
                                        break
                                break
                    if snoRNA not in linked_snoRNA.keys():
                        # all transcripts scanned, but no introns overlap the sno
                        # Select the first transcript (if it exists) that contains the sno but that overlaps with an exon
                        if len(potential_trans) > 0:
                            linked_snoRNA[snoRNA] = potential_trans[0]
                            exons = exon_df[exon_df['transcript_name'] == potential_trans[0]]
                            intronic_sno.loc[intronic_sno['gene_id'] == snoRNA, 'gene_id_HG'] = potential_trans[0] + '_mod'
                            if potential_trans[0] not in exon_hg.keys():
                                exon_hg[potential_trans[0] + '_mod'] = exons
                        else:
                            # although host gene overlaps the snoRNA, no individual transcript overlaps entirely the snoRNA
                            # thus this sno should be considered as intergenic (e.g. some snoRNAs in Snhg14 in mouse)
                            false_intronic_sno.append(snoRNA)





        # Groupby snoRNAs by the same host gene (drop false intronic snoRNAs)
        intronic_sno = intronic_sno[~intronic_sno['gene_id'].isin(false_intronic_sno)]
        intronic_sno = intronic_sno.sort_values(by=['seqname', 'strand', 'start'])
        intronic_sno['biotype'] = 'snoRNA'
        cols_mod = ['seqname', 'start', 'end', 'gene_id', 'score2', 'strand', 'biotype']
        intronic_cluster_nb, exonic_cluster_nb = 1, 1  # where exonic is one transcript (1 big exon) containing multiple snoRNAs
        for i, group in intronic_sno.groupby('gene_id_HG'):
            hg_gene_id = pd.unique(group['gene_id_HG'])[0]
            # Select exons of chosen transcript of that HG
            select_exons = exon_hg[hg_gene_id].sort_values(by='start').reset_index(drop=True)
            select_exons['biotype'] = 'host_gene_exon'
            select_exons['score2'] = select_exons['score']

            """HG with one snoRNA"""
            if len(group) == 1:
                sno_id = list(group['gene_id'])[0]
                sno_start = list(group['start'])[0]
                sno_end = list(group['end'])[0]


                if len(select_exons) == 1:  # only one exon in that HG
                    intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_exonic'
                
                else:  # multiple exons in HG
                    select_exons[cols_mod].to_csv(f'{species}_temp_exon.bed', sep='\t', index=False, header=False)
                    group[cols_mod].to_csv(f'{species}_temp_sno_mono.bed', sep='\t', index=False, header=False)
                    exon_bed = BedTool(f'{species}_temp_exon.bed')
                    sno_bed_mono = BedTool(f'{species}_temp_sno_mono.bed')
                    intersect = sno_bed_mono.intersect(exon_bed, s=True, wa=True, wb=True).to_dataframe(
                                            names=cols_mod+[i+'_HG' for i in cols_mod])

                    if len(intersect) == 0: # True intronic snoRNAs
                        intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_intronic'
                
                    else: # sno overlapping with exon, find if other transcript might have non-overlapping exons with that sno
                        # Set as mono_exonic by default, might be changed below if another transcript it found
                        intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_exonic'
                        trans_id = pd.unique(select_exons['transcript_id'])[0]
                        select_hg = temp_hg_gtf[(temp_hg_gtf['gene_id'] == hg_gene_id) & 
                                    (temp_hg_gtf['transcript_id'] != trans_id) & (temp_hg_gtf['feature'] == 'exon')]
                        # Iterate through transcript and each of their exon pairs to see if a pair of exon might surround the snoRNA
                        for i, trans in select_hg.groupby('transcript_id'):
                            trans = trans.sort_values(by='start').reset_index(drop=True)
                            exons_pairs = [[trans.iloc[i], trans.iloc[i+1]] for i in range(len(trans) - 1)]
                            for pair in exons_pairs:
                                exon1, exon2 = pair[0], pair[1]
                                # If an intron is found overlapping that sno, change its location to intronic instead of exonic
                                if (exon1['end'] < sno_start) & (exon2['start'] > sno_end):
                                    intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_intronic'
                                    break
                            break
                    sp.call(f'rm {species}_temp_exon.bed {species}_temp_sno_mono.bed', shell=True)


            """HG with multiple snoRNAs"""
            if len(group) > 1:
                if len(select_exons) == 1: # For transcript (1 big exon) containing snoRNAs
                    temp_cluster = []
                    for sno_id in group['gene_id']:
                        sno_row = group[group['gene_id'] == sno_id]
                        if (list(sno_row['start'])[0] >= list(select_exons['start'])[0]) & (list(sno_row['end'])[0] <= list(select_exons['end'])[0]):
                            temp_cluster.append(sno_id)
                    intronic_sno.loc[intronic_sno['gene_id'].isin(temp_cluster), 'cluster'] = f'exonic_cluster_{exonic_cluster_nb}'
                    exonic_cluster_nb += 1
                else:  # multiple exons in HG
                    exon_pairs = [[select_exons.iloc[i], select_exons.iloc[i+1]] for i in range(len(select_exons) - 1)]
                    temp_intronic_cluster = {}
                    # Iterate through snoRNAs within the same HG
                    for sno_id in group['gene_id']:
                        sno_row = group[group['gene_id'] == sno_id]
                        for j, pair in enumerate(exon_pairs):
                            exon1, exon2 = pair[0], pair[1]
                            # Select only truly intronic snoRNAs (non-overlapping with exons)
                            if (list(sno_row['start'])[0] > exon1['end']) & (list(sno_row['end'])[0] < exon2['start']):
                                temp_intronic_cluster[sno_id] = j
                                break

                    # Groupby truly intronic snoRNAs by the intron in which they are (and might share with other sno)
                    grouped_dict = {}
                    for k, v in temp_intronic_cluster.items():
                        if v not in grouped_dict:
                            grouped_dict[v] = [k]
                        else:
                            grouped_dict[v].append(k)

                    # Select firstly multi-sno introns (intronic clusters)
                    multi_sno_introns = {k:v for k,v in grouped_dict.items() if len(v) > 1}
                    for k,v in multi_sno_introns.items():
                        intronic_sno.loc[intronic_sno['gene_id'].isin(v), 'cluster'] = f'intronic_cluster_{intronic_cluster_nb}'
                        intronic_cluster_nb += 1

                    # Select secondly mono-sno introns
                    mono_sno_introns = {k:v for k,v in grouped_dict.items() if len(v) == 1}
                    for k,v in mono_sno_introns.items():
                        intronic_sno.loc[intronic_sno['gene_id'].isin(v), 'cluster'] = 'mono_intronic'
                    
                    
                    # For sno potentially overlapping with exons
                    for sno_id in group['gene_id']:
                        sno_start = list(group[group['gene_id'] == sno_id]['start'])[0]
                        sno_end = list(group[group['gene_id'] == sno_id]['end'])[0]
                        if sno_id not in temp_intronic_cluster.keys():  # potentially exonic snoRNAs
                            # Find if other transcript might have non-overlapping exons with that sno
                            # Set as mono_exonic by default, might be changed below if another transcript it found
                            intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_exonic'
                            trans_id = pd.unique(select_exons['transcript_id'])[0]
                            select_hg = temp_hg_gtf[(temp_hg_gtf['gene_id'] == hg_gene_id) & 
                                        (temp_hg_gtf['transcript_id'] != trans_id) & (temp_hg_gtf['feature'] == 'exon')]
                            # Iterate through transcript and each of their exon pairs to see if a pair of exon might surround the snoRNA
                            for i, trans in select_hg.groupby('transcript_id'):
                                trans = trans.sort_values(by='start').reset_index(drop=True)
                                exons_pairs = [[trans.iloc[i], trans.iloc[i+1]] for i in range(len(trans) - 1)]
                                for pair in exons_pairs:
                                    exon1, exon2 = pair[0], pair[1]
                                    # If an intron is found overlapping that sno, change its location to intronic instead of exonic
                                    if (exon1['end'] < sno_start) & (exon2['start'] > sno_end):
                                        intronic_sno.loc[intronic_sno['gene_id'] == sno_id, 'cluster'] = 'mono_intronic'
                                        break
                                break



        # Define the localization of intronic snoRNAs
        intronic_sno['cluster'] = intronic_sno['cluster'].fillna('NOT_DEFINED_intronic')
        intronic_sno.loc[intronic_sno['cluster'].str.contains('intronic_cluster_'), 'localization'] = 'intronic_cluster'
        intronic_sno.loc[intronic_sno['cluster'].str.contains('exonic_cluster_'), 'localization'] = 'exonic_cluster'
        intronic_sno.loc[intronic_sno['cluster'] == 'mono_intronic', 'localization'] = 'mono_intronic'
        intronic_sno.loc[intronic_sno['cluster'] == 'mono_exonic', 'localization'] = 'mono_exonic'
        intronic_sno.loc[intronic_sno['cluster'] == 'NOT_DEFINED_intronic', 'localization'] = 'NOT_DEFINED_intronic'

        # Define the localization of fake intronic snoRNAs (included in a host gene, but not in any individual transcript of that HG)
        fake_intronic_sno = initial_intronic_sno[initial_intronic_sno['gene_id'].isin(false_intronic_sno)]
        fake_intronic_sno['cluster'] = 'intergenic'
        fake_intronic_sno['localization'] = 'mono_intergenic'
        fake_intronic_sno['gene_id_HG'] = None


        # Concat dfs
        if 'gene_id_HG' not in intergenic_sno.columns:
            intergenic_sno['gene_id_HG'] = None
        final_df = pd.concat([intergenic_sno[['gene_id', 'gene_id_HG', 'cluster', 'localization']], 
                    fake_intronic_sno[['gene_id', 'gene_id_HG', 'cluster', 'localization']],
                    intronic_sno[['gene_id', 'gene_id_HG', 'cluster', 'localization']]])
        print(coll.Counter(final_df['localization']))

        final_df.to_csv(snakemake.output.sno_cluster, sep='\t', index=False)

        sp.call(f'rm hg_gtf_temp_{species}.gtf', shell=True)

    else:  # only intergenic sno annotated in that species (i.e. Tetrahymena thermophila)
        intergenic_sno['gene_id_HG'] = None
        print(coll.Counter(intergenic_sno['localization']))

        intergenic_sno[['gene_id', 'gene_id_HG', 'cluster', 'localization']].to_csv(snakemake.output.sno_cluster, 
                                                                                sep='\t', index=False)
