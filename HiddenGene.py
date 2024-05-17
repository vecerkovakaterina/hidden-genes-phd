from collections import Counter

import polars as pl

from hidden_genes_phd.Genome import Genome
from hidden_genes_phd.OrthologyGroup import OrthologyGroup
from hidden_genes_phd.OrthologyTable import OrthologyTable


class HiddenGene:
    # hidden_genes_list = []

    def __init__(self, orthology_group, genome):
        self.orthology_group = orthology_group
        self.missing_from_genome = genome
        self.left_neighbor = None
        self.left_neighbor_orthology_group = None
        self.right_neighbor = None
        self.right_neighbor_orthology_group = None
        self.region_between_neighbors = None
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

    def get_neighbor_coordinates(self):
        left_neighbor_coordinates = (
            self.missing_from_genome.annotation.get_coordinates_by_ensembl_id(
                self.left_neighbor
            )
        )
        right_neighbor_coordinates = (
            self.missing_from_genome.annotation.get_coordinates_by_ensembl_id(
                self.right_neighbor
            )
        )
        return left_neighbor_coordinates, right_neighbor_coordinates

    @staticmethod
    def combine_left_right_region_sequence(left, right):
        left = left.split("\n")
        left_header = left.pop(0)
        left_sequence = "".join(left)

        right = right.split("\n")
        right_header = right.pop(0)
        right_sequence = "".join(right)

        combined_header = left_header + " + " + right_header
        combined_sequence = left_sequence + "N" + right_sequence
        combined = combined_header + "\n" + combined_sequence
        return combined

    def get_region_between_neighbors(self):
        (
            left_neighbor_coordinates,
            right_neighbor_coordinates,
        ) = self.get_neighbor_coordinates()

        left_neighbor_scaffold = Genome.access_ensembl_api(
            f"/lookup/id/{self.left_neighbor}?expand=1"
        )["seq_region_name"]
        right_neighbor_scaffold = Genome.access_ensembl_api(
            f"/lookup/id/{self.right_neighbor}?expand=1"
        )["seq_region_name"]
        scaffold = None
        if left_neighbor_scaffold == right_neighbor_scaffold:
            scaffold = left_neighbor_scaffold
        if scaffold and (left_neighbor_coordinates > right_neighbor_coordinates):
            tmp = left_neighbor_coordinates
            left_neighbor_coordinates = right_neighbor_coordinates
            right_neighbor_coordinates = tmp

        region_start = left_neighbor_coordinates[1]
        region_end = right_neighbor_coordinates[0]

        if scaffold:
            region_sequence_response = Genome.access_ensembl_api(
                f"/sequence/region/{self.missing_from_genome.annotation_name}/{scaffold}:{region_start}..{region_end}?",
                content_type="text/x-fasta",
            )
            if region_sequence_response.text:
                self.region_between_neighbors = region_sequence_response.text

        else:
            left_part_response = Genome.access_ensembl_api(
                f"/sequence/region/{self.missing_from_genome.annotation_name}/{left_neighbor_scaffold}:{region_start}..{region_start+10000000}?",
                content_type="text/x-fasta",
            )  # try to get maximum sequence to reach the end of the scaffold
            if left_part_response.status_code == 400:
                # TODO retrieve in pieces
                print("Trying to retrieve region longer than 10 Mbp!")
            right_part_response = Genome.access_ensembl_api(
                f"/sequence/region/{self.missing_from_genome.annotation_name}/{right_neighbor_scaffold}:{1}..{region_end}?",
                content_type="text/x-fasta",
            )
            if right_part_response.status_code == 400:
                # TODO retrieve in pieces
                print("Trying to retrieve region longer than 10 Mbp!")

            left_part = None
            right_part = None

            if left_part_response.text.startswith(">"):
                left_part = left_part_response.text
            if right_part_response.text.startswith(">"):
                right_part = right_part_response.text

            if left_part and right_part:
                region_sequence = self.combine_left_right_region_sequence(
                    left_part, right_part
                )
                self.region_between_neighbors = region_sequence

    def write_region_sequence_to_fasta(self):
        if self.region_between_neighbors:
            region_sequence_fasta_filename = Path(
                "regions_to_search",
                "_".join(
                    [
                        self.missing_from_genome.species_name,
                        self.left_neighbor,
                        self.right_neighbor,
                    ]
                ),
            )
            if not region_sequence_fasta_filename.is_file():
                with open(region_sequence_fasta_filename, "w") as fasta_file:
                    fasta_file.write(self.region_between_neighbors)
                self.region_between_neighbors_fasta = region_sequence_fasta_filename

    def find_overlapping_annotation(self):
        # TODO
        pass
