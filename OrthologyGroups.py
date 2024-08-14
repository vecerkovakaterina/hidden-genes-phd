class OrthologyGroups:
    def __init__(self):
        self.orthology_groups_list = []

    def __iter__(self):
        for group in self.orthology_groups_list:
            yield group

    def score_all_groups(self, orthology_table): # TODO include taxon score
        for group in self:
            group.assign_score(orthology_table)
        return self

    def remove_full_groups(self, orthology_table):
        for group in self:
            if group.is_full(orthology_table):
                self.orthology_groups_list.remove(group)
        return self

    def rank_groups_by_score(self, orthology_table):
        self.score_all_groups(orthology_table)
        self.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)
        return self

