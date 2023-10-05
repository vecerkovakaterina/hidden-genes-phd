import polars as pl


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


def get_max_possible_score(orthology_table, taxon=None):
    full_orthology_group = 1.0
    if taxon:
        full_orthology_group += (
            orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
            .select(pl.count())
            .item()
        )
    return full_orthology_group


def get_max_size_orthology_group(orthology_table):
    return len(orthology_table.orthology_taxonomy_df.columns) - 1


class OrthologyGroup:
    orthology_groups_list = []

    def __init__(self, orthologs):
        self.orthologs = orthologs
        self.score = None
        OrthologyGroup.orthology_groups_list.append(self)

    def calculate_score(self, orthology_table, taxon=None):
        max_orthologs = get_max_size_orthology_group(orthology_table)

        for i, og in enumerate(self.orthology_groups_list):
            number_orthologs = len([x for x in og.orthologs if str(x) != "nan"])
            og.score = number_orthologs / max_orthologs
            if taxon:
                og.score += taxon_score(i, orthology_table, taxon)

    def delete_full_groups(self, orthology_table):
        max_orthologs = get_max_size_orthology_group(orthology_table)
        for og in self.orthology_groups_list:
            if len(og.orthologs) == max_orthologs:
                self.orthology_groups_list.remove(og)

    @classmethod
    def rank_groups_by_score(cls):
        cls.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)
