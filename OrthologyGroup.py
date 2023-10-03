def taxon_score(col, ortho_table, taxon):
    species_in_taxon_df = ortho_table.orthology_taxonomy_df.loc[
        ortho_table.orthology_taxonomy_df["class"] == taxon
    ].copy()
    return species_in_taxon_df[col].notna().sum()


def get_max_possible_score(ortho_table_taxon, taxon=None):
    full_orthology_group = 1.0
    if taxon:
        full_orthology_group += len(
            ortho_table_taxon.loc[ortho_table_taxon["class"] == taxon].index
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
            og.score = len(og.orthologs) / max_orthologs
            if taxon:
                og.score += taxon_score(i, ortho_table, taxon)

    def delete_full_groups(self, orthology_table):
        max_orthologs = get_max_size_orthology_group(orthology_table)
        for og in self.orthology_groups_list:
            if len(og.orthologs) == max_orthologs:
                self.orthology_groups_list.remove(og)

    @classmethod
    def rank_groups_by_score(cls):
        cls.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)
