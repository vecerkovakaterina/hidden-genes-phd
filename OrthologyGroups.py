import concurrent.futures


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

    def remove_full_groups(self):
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(group.is_full): group for group in self}

            for future in as_completed(futures):
                group = futures[future]
                if future.result():
                    self.orthology_groups_list.remove(group)

        return self

    def rank_groups_by_score(self, orthology_table):
        self.score_all_groups(orthology_table)
        self.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)
        return self

    def create_custom_blast_dbs(self):
        for group in self:
            group.fasta_to_custom_blast_db()

    @staticmethod
    def write_ortholog_sequences_to_fasta(orthology_group):
        orthology_group.write_ortholog_sequences_to_fasta()

    def write_ortholog_sequences_to_fastas(self):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(OrthologyGroups.write_ortholog_sequences_to_fasta, self)

