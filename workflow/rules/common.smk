""" Common functions used across the workflow"""

def get_species_gtf(species):
    # Get the gtf of the genome of a given species
    species = str(species)
    protists = ["dictyostelium_discoideum", "giardia_lamblia", "tetrahymena_thermophila"]
    animals = ['ornithorhynchus_anatinus', 'gallus_gallus', 
                'caenorhabditis_elegans', 'mus_musculus']
    plants = ['arabidopsis_thaliana', 'oryza_sativa', 'triticum_aestivum']
    fungi = ['saccharomyces_cerevisiae', 'schizosaccharomyces_pombe', 
            'neurospora_crassa', 'candida_albicans']
    others = ['rattus_norvegicus', 'macaca_mulatta', 'danio_rerio', 'xenopus_tropicalis']
    if species in fungi:
        path = rules.download_fungi_gtf.output.gtf
    elif species == 'homo_sapiens':
        path = rules.download_human_gtf.output.gtf
    elif species == 'drosophila_melanogaster':
        path = rules.download_droso_gtf.output.gtf
    elif species in others:
        path = rules.download_other_gtf.output.gtf
    elif species in protists:
        path = rules.download_protist_gtf.output.gtf
    elif species in plants:
        path = rules.download_plant_gtf.output.gtf
    elif species in animals:
        path = rules.download_animal_gtf.output.gtf
    return path

def join_list(l, subl, remove=False):
    # From a list l, return a string of all items in subl joined by '|'
    small_list = [a for a in l if a in subl]
    # If we want to remove (instead of only keeping) items of subl from l
    if remove==True:
        small_list = [a for a in l if a not in subl]
    return "{}".format("|".join(small_list))

