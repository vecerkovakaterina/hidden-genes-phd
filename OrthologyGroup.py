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


class OrthologyGroup:
    orthology_groups_list = []

    def __init__(self, orthologs):
        self.orthologs = orthologs
        self.score = None
        OrthologyGroup.orthology_groups_list.append(self)

    def score(self, ortho_table, taxon=None):
        max_orthologs = len(ortho_table.orthology_taxonomy_df.columns) - 1

        for i, og in enumerate(self.orthology_groups_list):
            og.score = len(og.orthologs) / max_orthologs
            if taxon:
                og.score += taxon_score(i, ortho_table, taxon)

    def delete_full(self):
        # TODO delete ogs with all orthologs
        pass

    def rank(self):
        # TODO return list of ogs ranked by score
        pass
