from collections import Counter

import polars as pl

from hidden_genes_phd.OrthologyGroup import OrthologyGroup


class HiddenGene:
    # hidden_genes_list = []

    def __init__(self, orthology_group, genome):
        self.orthology_group = orthology_group
        self.missing_from_genome = genome
        self.left_neighbor = None
        self.right_neighbor = None
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

    def get_n_neighbors_side(self, genomes, n_neighbors_to_search=10, side=None):
        self.validate_side(side=side)
        neighbors_dfs = []
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
            if side == "left":
                neighbors_df = genomes[species].annotation.df.slice(
                    offset=index_of_ortholog_in_annotation - n_neighbors_to_search,
                    length=n_neighbors_to_search,
                )
            elif side == "right":
                neighbors_df = genomes[species].annotation.df.slice(
                    offset=index_of_ortholog_in_annotation - n_neighbors_to_search,
                    length=n_neighbors_to_search,
                )
            neighbors_dfs.append(neighbors_df)
        return neighbors_dfs

    def find_most_common_neighbor_side(self, neighbors_df, orthology_groups, side=None):
        # lists of neighbors as ensembl ids
        self.validate_side(side=side)
        neighbors_ensembl_ids = []
        for df in neighbors_df:
            neighbors_ensembl_ids.append(
                df.select(pl.col("ensembl_id")).to_series().to_list()
            )

        # reverse list if side left
        if side == "left":
            neighbors_ensembl_ids = [x[::-1] for x in neighbors_ensembl_ids]

        # ensembl ids to orthology groups
        for i, neighbors in enumerate(neighbors_ensembl_ids):
            for j, neighbor in enumerate(neighbors):
                for orthology_group in orthology_groups:
                    if neighbor in orthology_group.orthologs:
                        neighbors_ensembl_ids[i][j] = orthology_group

        counter = Counter(neighbors_ensembl_ids[0])
        for i in neighbors_ensembl_ids[1:]:
            counter.update(i)

        occurences = counter.most_common()
        if len(occurences) > 0:
            return occurences[0][0]
        else:
            return None

    # def set_neighbor_genes(self):
    #     self.left_neighbor = self.find_neighbor_gene()
    #     self.right_neighbor = self.find_neighbor_gene()

    def find_overlapping_annotation(self):
        # TODO
        pass
