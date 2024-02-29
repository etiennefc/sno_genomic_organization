rule get_sno_HG_bed:
    """ Create a bed of all snoRNAs in gtf and a bed for all other 
        genes that are not mncRNAs. Also return if they overlap 
        (i.e. if snoRNAs are intergenic or intronic in protein-coding 
        or non-coding host genes)."""
    input:
        gtf = 'data/references/gtf/{species}.gtf',
        human_snobed_snoDB = 'data/references/snoDB_BED.bed',
        human_snotype_snoDB = 'data/references/snoDB_data.tsv'
    output:
        sno_bed = 'results/bed/sno_{species}.bed',
        HG_bed = 'results/bed/potential_HG_{species}.bed',
        sno_pos = 'results/localization/sno_localization_{species}.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_sno_HG_bed.py"

rule sno_clusters:
    """ Determine if snoRNAs are in a cluster of snoRNAs 
        (within the same intron or as intergenic cluster)."""
    input:
        sno_bed = rules.get_sno_HG_bed.output.sno_bed,
        sno_pos = rules.get_sno_HG_bed.output.sno_pos,
        HG_bed = rules.get_sno_HG_bed.output.HG_bed,
        gtf = get_species_gtf
    output:
        sno_cluster = 'results/localization/sno_cluster_{species}.tsv',
    params:
        intergenic_distance = 500
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/sno_clusters.py"

rule filter_rnacentral_id:
    """ From Ensembl ids, get corresponding RNAcentral ids of snoRNAs and their taxon id. Create
        a URL for each snoRNA to search for their snoRNA type in their description in RNAcentral."""
    input:
        ids = rules.download_rnacentral_ids.output.ids,
        sno_bed = rules.get_sno_HG_bed.output.sno_bed
    output:
        urls = 'data/references/RNA_central_search/RNA_central_URL_{species}.txt',
        id_table = 'data/references/RNA_central_search/RNA_central_to_external_id_{species}.tsv'
    params:
        url_part = 'https://rnacentral.org/rna/'
    wildcard_constraints:
        species=join_list(config['species'], ['homo_sapiens', 'saccharomyces_cerevisiae'], remove=True)
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_rnacentral_id.py"

rule get_missing_sno_type:
    """ Fetch from RNAcentral the snoRNA type in the snoRNA description page. 
        If not present, run infernal with Rfam covariance models to find a 
        potential snoRNA family (and snoRNA type). Do this for all species 
        except human and S. cerevisiae because they have their own specialized 
        snoRNA database which curates their snoRNA type."""
    input:
        urls = rules.filter_rnacentral_id.output.urls,
        id_table = rules.filter_rnacentral_id.output.id_table,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm,
        rfam_sno_type_family_id = 'data/references/snoRNA_type_rfam_families.tsv',
        sno_bed = rules.get_sno_HG_bed.output.sno_bed,
        genome_fa = 'data/references/genome_fa/{species}_genome.fa'
    output:
        sno_type = "data/references/RNA_central_search/sno_type_{species}.tsv"
    wildcard_constraints:
        species=join_list(config['species'], ['homo_sapiens', 'saccharomyces_cerevisiae'], remove=True)
    conda:
        "../envs/infernal.yaml"
    script:
        "../scripts/python/get_missing_sno_type.py"

rule regroup_sno_type_info:
    """ Get the sno_type info for all species into a df."""
    input:
        rnacentral_sno_type = expand(rules.get_missing_sno_type.output.sno_type, 
                                species=[s for s in config['species'] if s not in 
                                ['homo_sapiens', 'saccharomyces_cerevisiae']]),
        sno_type_human = 'data/references/snoDB_data.tsv',
        sno_type_cerevisiae = 'data/references/saccharomyces_cerevisiae_snotype_umass.tsv',
        sno_type_other = 'data/references/dieci_et_al_2009_Genomics_summary.tsv'
    output:
        sno_type_all = 'data/references/filtered_sno_type/all_species_snoRNA_type.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/regroup_sno_type_info.py"



