#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
from gtfparse import read_gtf
import collections as coll

cols = ['seqname', 'start', 'end', 'gene_id', 'score2', 'strand', 'gene_name']
species = str(snakemake.wildcards.species)
print(species)

def drop_duplicates_biotype(group):
            # Drop duplicate rows by keeping in priority one row with protein_coding as biotype_HG, 
            # otherwise keep the second row that is not protein_coding 
            # Check if 'protein_coding' is present in the 'gene_biotype' column of that group in groupby
            if 'protein_coding' in group['gene_biotype_HG'].values:
                # return one occurence (row) with protein_coding as biotype_HG
                return group.loc[group['gene_biotype_HG'].eq('protein_coding').idxmax()]
            else:
                # return one occurence (row) with another biotype HG than protein_coding
                return group.loc[group['gene_biotype_HG'].ne('protein_coding').idxmax()]

if species != 'homo_sapiens':
    if species == 'dictyostelium_discoideum':
       # Select snoRNA lines in the gtf 
        sp.call(f"""awk '$0 !~ /#/ && $3=="gene" && $0 ~ /sno/' {snakemake.input.gtf} > temp_sno_{species}.gtf""", shell=True)
        sno_gtf = read_gtf(f'temp_sno_{species}.gtf')
        sno_gtf['score2'] = '.'
    else:
        # Select snoRNA lines in the gtf
        sp.call(f"""awk '$0 !~ /#/ && $3=="gene" && $0 ~ /gene_biotype "snoRNA"/' {snakemake.input.gtf} > temp_sno_{species}.gtf""", shell=True)
        sno_gtf = read_gtf(f'temp_sno_{species}.gtf')
        sno_gtf['score2'] = '.'

    # Create bed of snoRNA
    if len(sno_gtf) > 0:
        if 'gene_name' not in sno_gtf.columns:
            sno_gtf['gene_name'] = sno_gtf['gene_id']
        sno_gtf['gene_name'] = sno_gtf['gene_name'].fillna(sno_gtf['gene_id'])
        sno_gtf.loc[sno_gtf['gene_name'] == '', 'gene_name'] = sno_gtf['gene_id']
        sno_bed_df = sno_gtf[cols]
        sno_bed_df = sno_bed_df.drop_duplicates(subset=['seqname', 'start', 'end', 'strand'])
        sno_bed_df.to_csv(snakemake.output.sno_bed, sep='\t', index=False, header=False)

        # Create bed of potential host gene of snoRNA by excluding all mid-size noncoding RNAs as HG
        if species == "dictyostelium_discoideum":
            # Don't count "ncRNA" because that's how snoRNAs are annotated in that species
            mncRNA = ['ncRNA', 'miRNA', 'Mt_tRNA', 'ribozyme', 'scaRNA', 'scRNA', 'snoRNA', 'snRNA', 'sRNA', 'tRNA', 'vault_RNA', 'TEC', 'artifact']
        else:
            mncRNA = ['miRNA', 'Mt_tRNA', 'ribozyme', 'scaRNA', 'scRNA', 'snoRNA', 'snRNA', 'sRNA', 'tRNA', 'vault_RNA', 'TEC', 'artifact']
        mncRNA = '|'.join([f'gene_biotype "{i}"' for i in mncRNA])
        sp.call(f"""awk '$0 !~ /#/ && $3=="gene" && $0 !~ /{mncRNA}/' {snakemake.input.gtf} > temp_HG_{species}.gtf""", shell=True)
        HG_gtf = read_gtf(f'temp_HG_{species}.gtf')
        HG_gtf['score2'] = '.'
        if 'gene_name' not in HG_gtf.columns:
            HG_gtf['gene_name'] = HG_gtf['gene_id']
        HG_bed = HG_gtf[cols + ['gene_biotype']]
        HG_bed.to_csv(snakemake.output.HG_bed, sep='\t', index=False, header=False)


        # Load the sno and HG beds
        sno_bed = BedTool(snakemake.output.sno_bed)
        HG_bed = BedTool(snakemake.output.HG_bed)

        # Get the intersect of snoRNA and HG to find the nb of intronic snoRNAs
        intronic_sno = sno_bed.intersect(HG_bed, s=True, f=1, wa=True, wb=True).to_dataframe(
                                        names=cols+[i+'_HG' for i in cols]+['gene_biotype_HG'])


        # Drop duplicate rows based on sno gene_id
        if len(intronic_sno) > 0:
            intronic_sno = intronic_sno.groupby('gene_id').apply(drop_duplicates_biotype).reset_index(drop=True)
            print(coll.Counter(intronic_sno.gene_biotype_HG))

            # Replace all HG biotypes that are not protein_coding to non_coding
            intronic_sno.loc[intronic_sno['gene_biotype_HG'] != 'protein_coding', 'gene_biotype_HG'] = 'non_coding'

            # Merge intronic_sno to intergenic sno
            all_sno_df = sno_bed_df.merge(intronic_sno[['gene_id', 'gene_id_HG', 'gene_name_HG', 'gene_biotype_HG']], 
                                how='left', on='gene_id')
            all_sno_df['gene_biotype_HG'] = all_sno_df['gene_biotype_HG'].fillna('intergenic')
            all_sno_df.to_csv(snakemake.output.sno_pos, index=False, sep='\t')
        else:
            all_sno_df = sno_bed_df.copy()
            all_sno_df['gene_biotype_HG'] = 'intergenic'
            all_sno_df.to_csv(snakemake.output.sno_pos, index=False, sep='\t')
        print(all_sno_df)

        sp.call(f'rm temp_sno_{species}.gtf temp_HG_{species}.gtf', shell=True)

    else:
        # No snoRNA annotated in that species gtf
        pd.DataFrame().to_csv(snakemake.output.sno_bed, sep='\t', index=False, header=False)
        pd.DataFrame().to_csv(snakemake.output.HG_bed, sep='\t', index=False, header=False)
        pd.DataFrame().to_csv(snakemake.output.sno_pos, index=False, sep='\t')
        sp.call(f'rm temp_sno_{species}.gtf', shell=True)
else:
    # Deal with snoRNAs in snoDB (already merged no duplicates)
    sno_bed = pd.read_csv(snakemake.input.human_snobed_snoDB, sep='\t', 
                names=['seqname', 'start', 'end', 'gene_id', 'score2', 'strand'])
    sno_bed['seqname'] = sno_bed['seqname'].str.strip('chr')
    sno_type = pd.read_csv(snakemake.input.human_snotype_snoDB, sep='\t')
    sno_type = sno_type.rename(columns={'snoDB ID': 'gene_id', 'Symbol': 'gene_name'})

    # Select only snoRNAs of type C/D, H/ACA, AluACA or unknown
    sno_ids = list(sno_type[sno_type['Box Type'].isin(['H/ACA', 'C/D', 'AluACA', 'unknown'])]['gene_id'])
    sno_bed = sno_bed[sno_bed['gene_id'].isin(sno_ids)]
    sno_bed = sno_bed.merge(sno_type[['gene_id', 'gene_name']], how='left', on='gene_id')
    sno_bed['gene_name'] = sno_bed['gene_name'].fillna(sno_bed['gene_id'])
    sno_bed_df = sno_bed.copy()
    sno_bed.to_csv(snakemake.output.sno_bed, sep='\t', index=False, header=False)

    # Create bed of potential host gene of snoRNA by excluding all mid-size noncoding RNAs as HG
    mncRNA = ['miRNA', 'Mt_tRNA', 'ribozyme', 'scaRNA', 'scRNA', 'snoRNA', 'snRNA', 'sRNA', 'tRNA', 'vault_RNA', 'TEC', 'artifact']
    mncRNA = '|'.join([f'gene_biotype "{i}"' for i in mncRNA])
    sp.call(f"""awk '$0 !~ /#/ && $3=="gene" && $0 !~ /{mncRNA}/' {snakemake.input.gtf} > temp_HG_{species}.gtf""", shell=True)
    HG_gtf = read_gtf(f'temp_HG_{species}.gtf')
    HG_gtf['score2'] = '.'
    if 'gene_name' not in HG_gtf.columns:
        HG_gtf['gene_name'] = HG_gtf['gene_id']
    HG_bed = HG_gtf[cols + ['gene_biotype']]
    HG_bed.to_csv(snakemake.output.HG_bed, sep='\t', index=False, header=False)


    # Load the sno and HG beds
    sno_bed = BedTool(snakemake.output.sno_bed)
    HG_bed = BedTool(snakemake.output.HG_bed)

    # Get the intersect of snoRNA and HG to find the nb of intronic snoRNAs
    intronic_sno = sno_bed.intersect(HG_bed, s=True, f=1, wa=True, wb=True).to_dataframe(
                                    names=cols+[i+'_HG' for i in cols]+['gene_biotype_HG'])


    # Drop duplicate rows based on sno gene_id
    if len(intronic_sno) > 0:
        intronic_sno = intronic_sno.groupby('gene_id').apply(drop_duplicates_biotype).reset_index(drop=True)
        print(coll.Counter(intronic_sno.gene_biotype_HG))

        # Replace all HG biotypes that are not protein_coding to non_coding
        intronic_sno.loc[intronic_sno['gene_biotype_HG'] != 'protein_coding', 'gene_biotype_HG'] = 'non_coding'

        # Merge intronic_sno to intergenic sno
        all_sno_df = sno_bed_df.merge(intronic_sno[['gene_id', 'gene_id_HG', 'gene_name_HG', 'gene_biotype_HG']], 
                            how='left', on='gene_id')
        all_sno_df['gene_biotype_HG'] = all_sno_df['gene_biotype_HG'].fillna('intergenic')
        all_sno_df.to_csv(snakemake.output.sno_pos, index=False, sep='\t')


    sp.call(f'rm temp_HG_{species}.gtf', shell=True)






