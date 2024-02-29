#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import collections as coll
from pybedtools import BedTool
import glob 
import re

species = str(snakemake.wildcards.species)
urls = snakemake.input.urls
id_table = pd.read_csv(snakemake.input.id_table, sep='\t')
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', 
                names=['chr', 'start', 'end', 'gene_id', 
                        'score', 'strand', 'gene_name'])
rfam_snoRNA_type = pd.read_csv(snakemake.input.rfam_sno_type_family_id, sep='\t')
cov_mod = snakemake.input.rfam_cm
genome = snakemake.input.genome_fa
cols = ['gene_id', 'gene_name', 'snoRNA_type']

if species == 'giardia_lamblia':
    pd.DataFrame(columns=cols).to_csv(snakemake.output.sno_type, sep='\t', index=False)

else:
    # Format gene_id column to have the right external gene_id
    source = pd.unique(id_table['source'])[0]
    taxon_id = pd.unique(id_table['taxon_id'])[0]
    if source in ['ENSEMBL_PLANTS', 'FLYBASE', 'ENSEMBL_FUNGI', 'POMBASE']:
        id_table['gene_id_mod'] = id_table['gene_id2']
    if source in ['DICTYBASE', 'WORMBASE', 'ENSEMBL_PROTISTS']:
        id_table['gene_id_mod'] = id_table['gene_id']
    if source == 'ENSEMBL':
        id_table['gene_id_mod'] = id_table['gene_id2'].replace(
                                    to_replace=r'()\.[0-9]*$', 
                                    value=r'\1', regex=True)

    id_table = id_table.merge(sno_bed[['chr', 'start', 'end', 'gene_id', 'strand']], how='left', left_on='gene_id_mod', right_on='gene_id')
    id_table = id_table.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])

    # Create dict of snoRNA id and their respective rnacentral id
    d = dict(zip(id_table.gene_id_mod, id_table.rnacentral_id))
    d = {k:f'{v}_{taxon_id}' for k,v in d.items()}
    
    # Fetch snoRNA description on RNAcentral
    sp.call(f'wget --quiet -i {urls}', shell=True)

    # Grep snoRNA type for each snoRNA
    sno_desc_files = glob.glob('URS*')
    cd_sno, haca_sno, unknown_sno = [], [], []
    unknown_sno_seq = {}
    rnacentral_id_no_seq = []
    for sno in d.keys():
        description, sequence = '', ''
        search_file = [p for p in sno_desc_files if d[sno] in p][0]
        with open(search_file, 'r') as f:
            for line in f:
                if 'description' in line:
                    description += line
                if 'rna-sequence' in line:
                    sequence += line
        
        cd_motif = 'SNORD|Snord|U3|U8|U13|C/D|C_D'
        haca_motif = 'SNORA|Snora|H/ACA|HACA|H_ACA'
        cd = re.findall(cd_motif, description)
        haca = re.findall(haca_motif, description)
        if (len(cd) > 0) & (len(haca) == 0):
            cd_sno.append(sno)
        elif (len(haca) > 0) & (len(cd) == 0):
            haca_sno.append(sno)
        elif (len(haca) > 0) & (len(cd) > 0):
            print(sno, d[sno], 'BOTHHHHHHH')
            unknown_sno.append(sno)
        else:
            #print(sno, d[sno], 'UNKNOWN')
            if sequence != '':
                seq = sequence.split('"rna-sequence">')[1].split('<')[0]
                unknown_sno_seq[sno] = seq
                unknown_sno.append(sno)
            else:
                rnacentral_id_no_seq.append(sno)
    rnacentral_id_no_seq_df = sno_bed[sno_bed['gene_id'].isin(rnacentral_id_no_seq)]
    sno_wo_rnacentral = sno_bed[~sno_bed['gene_id'].isin(id_table['gene_id_mod'])]
    sno_wo_rnacentral = sno_wo_rnacentral[~sno_wo_rnacentral['gene_id'].isin(unknown_sno)]
    print('\nAfter RNAcentral search:')
    print(f'CD: {len(cd_sno)}, HACA: {len(haca_sno)}, UNKNOWN:{len(sno_bed) - len(cd_sno) - len(haca_sno)}\n')
    
    # Get the sequence of the snoRNA without RNAcentral ids
    if (len(sno_wo_rnacentral) > 0) | (len(rnacentral_id_no_seq_df) > 0):
        if (len(sno_wo_rnacentral) > 0) & (len(rnacentral_id_no_seq_df) > 0):
            concat_df = pd.concat([sno_wo_rnacentral, rnacentral_id_no_seq_df])
            concat_df.to_csv(f'sno_wo_rnacentral_{species}.bed', sep='\t', index=False, header=False)
            sno_wo_rnacentral_bed = BedTool(f'sno_wo_rnacentral_{species}.bed')
            fasta = sno_wo_rnacentral_bed.sequence(fi=genome, nameOnly=True, s=True)
        if (len(sno_wo_rnacentral) > 0) & (len(rnacentral_id_no_seq_df) == 0):
            sno_wo_rnacentral.to_csv(f'sno_wo_rnacentral_{species}.bed', sep='\t', index=False, header=False)
            sno_wo_rnacentral_bed = BedTool(f'sno_wo_rnacentral_{species}.bed')
            fasta = sno_wo_rnacentral_bed.sequence(fi=genome, nameOnly=True, s=True)
        if (len(sno_wo_rnacentral) == 0) & (len(rnacentral_id_no_seq_df) > 0):
            rnacentral_id_no_seq_df.to_csv(f'sno_wo_rnacentral_{species}.bed', sep='\t', index=False, header=False)
            sno_wo_rnacentral_bed = BedTool(f'sno_wo_rnacentral_{species}.bed')
            fasta = sno_wo_rnacentral_bed.sequence(fi=genome, nameOnly=True, s=True)
        with open(fasta.seqfn, 'r') as fasta_file:
            for line in fasta_file:
                if '>' in line:
                    name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                else:
                    seq2 = line.strip('\n')
                    unknown_sno_seq[name] = seq2
        sp.call(f'rm sno_wo_rnacentral_{species}.bed', shell=True)

    # For the unknown snoRNA type, use infernal and Rfam covariance models to 
    # predict their family and thus their type
    with open(f'temp_unknown_sno_{species}.fa', 'w') as f:
        for sno, seq_ in unknown_sno_seq.items():
            f.write(f'>{sno}\n{seq_}\n') 
    
    sp.call(f'cmpress -F {cov_mod} && cmscan --cut_ga --rfam --nohmmonly -o infernal_{species}.aln --tblout infernal_{species}.tblout {cov_mod} temp_unknown_sno_{species}.fa', shell=True)

    # Filter infernal output to find the snoRNA type associated to a snoRNA family
    sp.call(f"grep -v '#' infernal_{species}.tblout | awk -v OFS='\t' '{{print $3, $2}}' > filtered_infernal_{species}.tsv", shell=True)
    filtered_infernal = pd.read_csv(f'filtered_infernal_{species}.tsv', sep='\t', names=['gene_id', 'Rfam_family_id'])
    filtered_infernal = filtered_infernal.merge(rfam_snoRNA_type, how='left', on='Rfam_family_id')
    filtered_infernal['snoRNA_type'] = filtered_infernal['snoRNA_type'].fillna('Unknown')
    filtered_infernal = filtered_infernal.drop_duplicates(subset='gene_id')
    print(f'\nAfter infernal Rfam models search, out of the {len(sno_bed) - len(cd_sno) - len(haca_sno)} UNKNOWN:')
    print(f'CD: {len(filtered_infernal[filtered_infernal["snoRNA_type"] == "C/D"])}, HACA: {len(filtered_infernal[filtered_infernal["snoRNA_type"] == "H/ACA"])}, UNKNOWN:{len(sno_bed) - len(cd_sno) - len(haca_sno) - len(filtered_infernal[filtered_infernal["snoRNA_type"].isin(["H/ACA", "C/D"])])}\n')

    
    # Get snoRNA type in a proper column
    final_df = sno_bed[['gene_id', 'gene_name']].merge(filtered_infernal, how='left', on='gene_id')
    final_df.loc[final_df['gene_id'].isin(cd_sno), 'snoRNA_type'] = 'C/D'
    final_df.loc[final_df['gene_id'].isin(haca_sno), 'snoRNA_type'] = 'H/ACA'
    final_df = final_df.drop(columns=['Rfam_family_id'])
    final_df['snoRNA_type'] = final_df['snoRNA_type'].fillna('Unknown')

    final_df.to_csv(snakemake.output.sno_type, sep='\t', index=False)
    print('Final snoRNA type number:', coll.Counter(final_df['snoRNA_type']))

    sp.call('rm URS*', shell=True)
    sp.call(f'rm temp_unknown_sno_{species}.fa infernal_{species}.aln infernal_{species}.tblout filtered_infernal_{species}.tsv', shell=True)

    
