from pathlib import Path


from OrthologyTable import OrthologyTable
from HiddenGene import HiddenGene
from Genomes import Genomes
from OrthologyGroups import OrthologyGroups


target_taxon = "Mammalia"
input_orthology_table = Path("ten_columns_orthology_table.csv")

# parse orthology table
orthology_table = OrthologyTable(input_orthology_table)

# create Genome instances
genomes = Genomes()
orthology_table.create_genomes_for_species_in_table(genomes)

# add taxonomy class to each species
orthology_table = orthology_table.add_taxonomy_class_to_df(genomes)

# create orthology groups
orthology_groups = OrthologyGroups()
orthology_table.create_orthology_groups(orthology_groups)

# create custom blast dbs for orthology groups
orthology_groups.write_ortholog_sequences_to_fastas()
orthology_groups.create_custom_blast_dbs()


# create hidden gene instances
hidden_genes = []
species_names = orthology_table.get_annotation_names()
for i, group in enumerate(orthology_groups):
    for j, ortholog in enumerate(group.orthologs):
        if ortholog == "nan":
            hidden_genes.append(
                HiddenGene(group, genomes.genomes_dict[species_names[j]])
            )

# find most common left and right neighbor
for hidden_gene in hidden_genes[:10]:
    side = "left"
    each_side_search = 5
    left_neighbors = hidden_gene.get_n_neighbors_side(
        genomes, orthology_table, n_neighbors_to_search=each_side_search, side=side
    )
    if not left_neighbors:
        continue
    hidden_gene.left_neighbor_orthology_group = (
        hidden_gene.find_most_common_neighbor_side(
            left_neighbors, orthology_groups, side=side
        )
    )

    side = "right"
    right_neighbors = hidden_gene.get_n_neighbors_side(
        genomes, orthology_table, n_neighbors_to_search=each_side_search, side=side
    )
    hidden_gene.right_neighbor_orthology_group = (
        hidden_gene.find_most_common_neighbor_side(
            right_neighbors, orthology_groups, side=side
        )
    )

# remove groups with all orthologs
orthology_groups = orthology_groups.remove_full_groups()

# rank the orthology groups by their score
orthology_groups = orthology_groups.rank_groups_by_score(orthology_table)


for hidden_gene in hidden_genes[:100]:
    print(hidden_gene.missing_from_genome.species_name)
    hidden_gene.assign_neighbor_ensembl_id(orthology_table)
    print(
        f"Hidden gene left neighbor group: {hidden_gene.left_neighbor_orthology_group}, right nehighbor group: {hidden_gene.right_neighbor_orthology_group}"
    )
    print(
        f"Hidden gene left neighbor: {hidden_gene.left_neighbor}, right nehighbor: {hidden_gene.right_neighbor}"
    )
    if hidden_gene.left_neighbor and hidden_gene.right_neighbor:
        try:
            hidden_gene.get_region_between_neighbors()
            hidden_gene.write_region_sequence_to_fasta()
            hidden_gene.blast_region_to_ortholog_database()
            # TODO if the region between neighbors successfully retrieved, blast it against the library
        except KeyError:
            print("Unknown scaffold of neighbor(s)!")
            continue
        # print(
        #     orthology_table.get_species_ortholog_ensmebl_id_from_group(
        #         hidden_gene.missing_from_genome.species_name,
        #         hidden_gene.left_neighbor_orthology_group,
        #     )
        # )
