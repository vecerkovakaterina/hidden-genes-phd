import csv

from HiddenGene import HiddenGene


class HiddenGenes:
    def __init__(self):
        self.hidden_genes_list = []

    def __iter__(self):
        for hidden_gene in self.hidden_genes_list:
            yield hidden_gene

    def add_hidden_gene(self, hidden_gene):
        self.hidden_genes_list.append(hidden_gene)

    def create_hidden_genes_list(self, orthology_table, orthology_groups, genomes):
        species_names = orthology_table.get_annotation_names()
        for i, group in enumerate(orthology_groups):
            for j, ortholog in enumerate(group.orthologs):
                if ortholog == "nan":
                    self.add_hidden_gene(
                        HiddenGene(group, genomes.genomes_dict[species_names[j]])
                    )

    def export_to_cvs(self, file_path):
        attributes = [
            "missing_from_genome_name",
            "orthology_group_first_ortholog",
            "left_neighbor",
            "right_neighbor",
            "region_between_neighbors",
            "region_between_neighbors_fasta",
            "blast_output",
            "coordinates",
            "sequence",
            "overlaps_with",
        ]

        with open(file_path, mode="w") as file:
            writer = csv.DictWriter(file, fieldnames=attributes)
            writer.writeheader()
            for hidden_gene in self.hidden_genes_list:
                writer.writerow(
                    {
                        "missing_from_genome_name": hidden_gene.missing_from_genome.species_name,
                        "orthology_group_first_ortholog": next(
                            iter(hidden_gene.orthology_group.orthologs), None
                        ),
                        "left_neighbor": hidden_gene.left_neighbor,
                        "right_neighbor": hidden_gene.right_neighbor,
                        "region_between_neighbors": hidden_gene.region_between_neighbors,
                        "region_between_neighbors_fasta": hidden_gene.region_between_neighbors_fasta,
                        "blast_output": hidden_gene.blast_output,
                        "coordinates": hidden_gene.coordinates,
                        "sequence": hidden_gene.sequence,
                        "overlaps_with": hidden_gene.overlaps_with,
                    }
                )
