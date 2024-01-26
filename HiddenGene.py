from collections import Counter

import polars as pl

from hidden_genes_phd.OrthologyGroup import OrthologyGroup


class HiddenGene:
    # hidden_genes_list = []

    def __init__(self, orthology_group, genome):
        self.orthology_group = orthology_group
        self.missing_from_genome = genome
        self.left_neighbor = None
        self.left_neighbor_orthology_group = None
        self.right_neighbor = None
        self.right_neighbor_orthology_group = None
        self.coordinates = None
        self.sequence = None
        self.overlaps_with = None

        # self.missing_from_annotation = None
        # self.missing_from_assembly = None
        # self.missing_from_sequence = None

        # HiddenGene.hidden_genes_list.append(self)

    @staticmethod
    def validate_side(side=None):
        valid_sides = ["left", "right"]
        if side not in valid_sides:
            raise ValueError(f"Side must be of {valid_sides}")

    def get_n_neighbors_side(
        self, genomes, orthology_table, n_neighbors_to_search=10, side=None
    ):
        self.validate_side(side=side)
        neighbors_lists = []
        for i, species in enumerate(self.orthology_group.species_with_ortholog):
            orthologs = OrthologyGroup.drop_nans_from_list(
                self.orthology_group.orthologs
            )
            find_ortholog_in_annotation = pl.col("ensembl_id") == orthologs[i]
            index_of_ortholog_in_annotation = (
                genomes[species]
                .annotation.df.filter(find_ortholog_in_annotation)
                .select(pl.col("index"))
                .item()
            )

            neighbors = []
            counter = 1
            while len(neighbors) < n_neighbors_to_search:
                if side == "left":
                    next_gene_ensembl_id = (
                        genomes[species]
                        .annotation.df.slice(
                            offset=index_of_ortholog_in_annotation - counter,
                            length=1,
                        )
                        .select(pl.col("ensembl_id"))
                        .item()
                    )
                elif side == "right":
                    next_gene_ensembl_id = (
                        genomes[species]
                        .annotation.df.slice(
                            offset=index_of_ortholog_in_annotation + counter,
                            length=1,
                        )
                        .select(pl.col("ensembl_id"))
                        .item()
                    )
                # check if next gene present in annotation with hidden gene
                if (
                    orthology_table.get_ensembl_id_of_ortholog_in_species(
                        next_gene_ensembl_id,
                        species,
                        self.missing_from_genome.species_name,
                    )
                    != "nan"
                ):
                    neighbors.append(next_gene_ensembl_id)
                counter += 1
            neighbors_lists.append(neighbors)

        return neighbors_lists

    def find_most_common_neighbor_side(
        self, neighbors_ensembl_ids, orthology_groups, side=None
    ):
        # reverse list if side left
        if side == "left":
            neighbors_ensembl_ids = [x[::-1] for x in neighbors_ensembl_ids]

        # ensembl ids to orthology groups
        for i, neighbors in enumerate(neighbors_ensembl_ids):
            for j, neighbor in enumerate(neighbors):
                for orthology_group in orthology_groups:
                    if neighbor in orthology_group.orthologs:
                        neighbors_ensembl_ids[i][j] = orthology_group
                if not isinstance(
                    neighbors_ensembl_ids[i][j],
                    OrthologyGroup,
                ):
                    neighbors_ensembl_ids[i][j] = None

        counter = Counter(neighbors_ensembl_ids[0])
        for i in neighbors_ensembl_ids[1:]:
            counter.update(i)

        occurences = (
            counter.most_common()
        )  # TODO more equally frequent ones, can check for the closest one from the most frequent ones
        if len(occurences) > 0:
            # take another one if the most frequent one is None
            if occurences[0][0] is None and len(occurences) > 1:
                return occurences[1][0]
            else:
                return occurences[0][0]
        else:
            return None

    def assign_neighbor_ensembl_id(self, orthology_table):
        self.left_neighbor = orthology_table.get_species_ortholog_ensmebl_id_from_group(
            self.missing_from_genome,
            self.left_neighbor_orthology_group,
        )
        self.right_neighbor = (
            orthology_table.get_species_ortholog_ensmebl_id_from_group(
                self.missing_from_genome,
                self.right_neighbor_orthology_group,
            )
        )

    def get_coordinates_of_neighbor(self):
        # TODO
        pass

    def get_region_between_neighbors(self):
        # TODO
        pass

    def find_overlapping_annotation(self):
        # TODO
        pass
