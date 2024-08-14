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

    @staticmethod
    def get_max_size_orthology_group(orthology_table):
        return len(orthology_table.orthology_taxonomy_df.columns) - 1

    def assign_score(self, orthology_table):
        max_orthologs = self.get_max_size_orthology_group(orthology_table)
        number_orthologs = len(self.drop_nans_from_orthologs_list())
        self.score = number_orthologs / max_orthologs

    def is_full(self, ortology_table):
        max_orthologs = self.get_max_size_orthology_group(ortology_table)
        if len(self.orthologs) == max_orthologs:
            return True
        else:
            return False

    def create_ortholog_sequences_list(self):
        orthologs = self.drop_nans_from_list(self.orthologs)
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
        aa_sequences_filename = Path(
            "ortholog_sequences",
            f"{'_'.join(self.drop_nans_from_list(self.orthologs)) + '.fa'}",
        )
        if not aa_sequences_filename.is_file():
            # check if file exists already
            aa_sequences = self.create_ortholog_sequences_list()
            aa_sequences = "\n".join(["\n".join(sequence) for sequence in aa_sequences])

            with open(aa_sequences_filename, "w") as fasta_file:
                fasta_file.write(aa_sequences)

            self.sequences_fasta = aa_sequences_filename

    def fasta_to_custom_blast_db(self):
        command = f"makeblastdb -in {self.sequences_fasta} -parse_seqids -dbtype prot"
        result = subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
        if result.returncode != 0:
            raise Exception(
                f"Error when creating custom BLAST db: {self.sequences_fasta}!"
            )
