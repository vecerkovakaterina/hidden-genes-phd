import subprocess
from pathlib import Path

import gget
import polars as pl


class OrthologyGroup:
    def __init__(self, orthologs, species_with_ortholog, orthology_groups):
        self.orthologs = orthologs
        self.species_with_ortholog = species_with_ortholog
        self.score = None
        self.sequences_fasta = None
        orthology_groups.orthology_groups_list.append(self)

    def __hash__(self):
        return hash(tuple(self.orthologs))

    def __eq__(self, other):
        if isinstance(other, OrthologyGroup):
            return self.orthologs == other.orthologs
        return False

    def drop_nans_from_orthologs_list(self):
        return [ortholog for ortholog in self.orthologs if str(ortholog) != "nan"]

    @staticmethod
    def taxon_score(col, orthology_table, taxon):
        nans_in_taxon = (
            orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
            .select(str(orthology_table.orthology_taxonomy_df.columns[col]))
            .null_count()
            .item()
        )
        species_in_taxon = (
            orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
            .select(str(orthology_table.orthology_taxonomy_df.columns[col]))
            .select(pl.count())
            .item()
        )
        orthologs_in_taxon = species_in_taxon - nans_in_taxon
        return orthologs_in_taxon

    @staticmethod
    def get_max_possible_score(orthology_table, taxon=None):
        full_orthology_group = 1.0
        if taxon:
            full_orthology_group += (
                orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
                .select(pl.count())
                .item()
            )
        return full_orthology_group

    def assign_score(self, orthology_table):
        max_orthologs = orthology_table.max_number_orthologs
        number_orthologs = len(self.drop_nans_from_orthologs_list())
        self.score = number_orthologs / max_orthologs
        return self

    def is_full(self):
        if "nan" in self.orthologs:
            return False
        else:
            return True

    def create_ortholog_sequences_list(self):
        orthologs = self.drop_nans_from_orthologs_list()
        ortholog_aa_sequences = []
        for ortholog in orthologs:
            try:
                sequence = gget.seq(ortholog, translate=True, verbose=False)
                if len(sequence) > 0:
                    ortholog_aa_sequences.append(sequence)
            except TypeError:
                continue

        return ortholog_aa_sequences

    def write_ortholog_sequences_to_fasta(self):
        orthologs = self.drop_nans_from_orthologs_list()
        if len(orthologs) > 0:
            aa_sequences_filename = Path(
                "ortholog_sequences",
                f"{'_'.join(orthologs) + '.fa'}",
            )
            if not aa_sequences_filename.is_file():
                aa_sequences = self.create_ortholog_sequences_list()

                if len(aa_sequences) > 0:
                    aa_sequences = "\n".join(
                        ["\n".join(sequence) for sequence in aa_sequences]
                    )

                    with open(aa_sequences_filename, "w") as fasta_file:
                        fasta_file.write(aa_sequences)

                        self.sequences_fasta = aa_sequences_filename

            else:
                self.sequences_fasta = aa_sequences_filename

    def fasta_to_custom_blast_db(self):
        if self.sequences_fasta is not None:
            command = (
                f"makeblastdb -in {self.sequences_fasta} -parse_seqids -dbtype prot"
            )
            result = subprocess.run(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode != 0:
                raise Exception(
                    f"Error when creating custom BLAST db: {self.sequences_fasta}!"
                )
