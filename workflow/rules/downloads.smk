rule download_human_gtf:
    """ Download gtf of human genome from Zenodo. Remove trailing tabs."""
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'  # Ensembl v110 supplemented with snoDB v2 and GtRNAdb
    params:
        link = config['download']['human_gtf']
    shell:
        "wget -O temp_gtf_homo_sapiens.gz {params.link} && "
        "gunzip temp_gtf_homo_sapiens.gz && "
        "mv temp_gtf_homo_sapiens {output.gtf} && "
        "sed -i 's/[\t]$//g' {output.gtf}"

rule download_sno_type_info:
    """ Download snoDB data (bed file and snoRNA type of human snoRNAs), 
        snoRNA type of S. cerevisiae (from Umass db), snoRNA type info 
        from the Dieci et al 2009 Genomics paper and snoRNA type per 
        Rfam family. The snotype per Rfam family was obtained manually 
        by saving (right-click) the Rfam webpage displaying all C/D 
        (and then H/ACA) snoRNA families and parsing the rfam family id 
        only using the following command: 
        grep ">RF" ~/Downloads/HACA_RFAM_families_2024.html | sed -E 's/.*">//g; s/<.*//g'"""
    output:
        snodb_bed = 'data/references/snoDB_BED.bed',
        sno_type_snodb = 'data/references/snoDB_data.tsv',
        sno_type_cerevisiae = 'data/references/saccharomyces_cerevisiae_snotype_umass.tsv', 
        sno_type_dieci = 'data/references/dieci_et_al_2009_Genomics_summary.tsv',
        sno_type_rfam = 'data/references/snoRNA_type_rfam_families.tsv'
    params:
        link1 = config['download']['snoDB_bed'],
        link2 = config['download']['snoDB_data'],
        link3 = config['download']['sno_type_cerevisiae'],
        link4 = config['download']['sno_type_dieci'],
        link5 = config['download']['sno_type_rfam']
    shell:
        "wget -O {output.snodb_bed} {params.link1} && "
        "wget -O {output.sno_type_snodb} {params.link2} && "
        "wget -O {output.sno_type_cerevisiae} {params.link3} && "
        "wget -O {output.sno_type_dieci} {params.link4} && "
        "wget -O {output.sno_type_rfam} {params.link5}"

rule download_rnacentral_ids:
    """ Get a table of conversion between RNA central ids and other expert 
        databases ids (ex: Ensembl, Wormbase, Pombase, Flybase, etc.)."""
    output:
        ids = 'data/references/rnacentral_ids.tsv'
    params:
        link = config['download']['rnacentral_id_conversion']
    shell:
        "wget -O temp_rnacentral.gz {params.link} && "
        "gunzip temp_rnacentral.gz && "
        "mv temp_rnacentral {output.ids}"

rule download_animal_gtf:
    """ Download the annotation (gtf) of different animals
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ['ornithorhynchus_anatinus', 
                    'gallus_gallus', 'caenorhabditis_elegans', 
                    'mus_musculus'])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_droso_gtf:
    """ Download the annotation (gtf) of drosophila
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ['drosophila_melanogaster'])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/gtf/{species}/*46.110.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_other_gtf:
    """ Download the annotation (gtf) of rat, macaque, zebrafish and frog
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ['rattus_norvegicus', 'macaca_mulatta', 'danio_rerio', 'xenopus_tropicalis'])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/gtf/{species}/*.110.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"


rule download_plant_gtf:
    """ Download the annotation (gtf) of different plants
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["oryza_sativa", "arabidopsis_thaliana", "triticum_aestivum"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_fungi_gtf:
    """ Download the reference genome (fasta file) of S. cerevisiae, 
        S. pombe, A. fumigatus, N. crassa and C. albicans
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae", 
                "schizosaccharomyces_pombe", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-58/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_protist_gtf:
    """ Download the annotation (gtf) of different protists 
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["dictyostelium_discoideum", 
                        "giardia_lamblia", "tetrahymena_thermophila"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-58/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_rfam_covariance_models:
    """ Download the Rfam library of covariance models that 
        will be used by Infernal."""
    output:
        rfam_cm = 'data/references/RFam.cm'
    params:
        link = config['download']['rfam_cm']
    shell:
        "wget {params.link} && "
        "gunzip Rfam.cm.gz && mv Rfam.cm {output.rfam_cm}"

rule download_mouse_genome:
    """ Download the reference genome (fasta file) of mouse
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["mus_musculus"])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/fasta/{species}/dna/*dna.primary_assembly.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_animal_genome:
    """ Download the reference genome (fasta file) of animals
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["macaca_mulatta", "rattus_norvegicus", "ornithorhynchus_anatinus",
                    "gallus_gallus", "danio_rerio", "xenopus_tropicalis", "caenorhabditis_elegans", 
                    "drosophila_melanogaster"])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_plant_genome:
    """ Download the reference genome (fasta file) of animals
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["arabidopsis_thaliana", "oryza_sativa", "triticum_aestivum"])
    params:
        link = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_protist_genome:
    """ Download the reference genome (fasta file) of protists species from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["tetrahymena_thermophila", 
              "dictyostelium_discoideum", "giardia_lamblia"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-58/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_yeast_genome:
    """ Download the reference genome (fasta file) of S. pombe, 
        N. crassa and C. albicans from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["schizosaccharomyces_pombe", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

