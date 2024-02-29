rule bar_sno_type_localization:
    """ Create a grouped stacked bar chart of the box type 
        and localization of snoRNAs per eukaryotic kingdoms."""
    input:
        sno_type = rules.regroup_sno_type_info.output.sno_type_all,
        sno_localization = expand(rules.sno_clusters.output.sno_cluster, **config),
        sno_pos = expand(rules.get_sno_HG_bed.output.sno_pos, **config)
    output:
        bar_animals = 'results/figures/sno_type_localization_animals.svg',
        bar_others = 'results/figures/sno_type_localization_others.svg',
        final_df = 'results/localization/final_table.tsv'
    params:
        colors_location_type = config['colors']['location_type']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/bar_sno_type_localization.py"

rule donut_host_biotype:
    """ Create a donut chart showing the distribution of host 
        gene biotype per species across eukaryotic kingdoms."""
    input:
        sno_localization = expand(rules.get_sno_HG_bed.output.sno_pos, **config)
    output:
        donut_animals = 'results/figures/donut_host_biotype_animals.svg',
        donut_others = 'results/figures/donut_host_biotype_others.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/donut_host_biotype.py"