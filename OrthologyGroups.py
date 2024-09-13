from concurrent.futures import ThreadPoolExecutor, as_completed


class OrthologyGroups:
    def __init__(self):
        self.orthology_groups_list = []

    def __iter__(self):
        for group in self.orthology_groups_list:
            yield group

    # def score_all_groups(self, orthology_table):  # TODO include taxon score
    #     for group in self:
    #         group.assign_score(orthology_table)
    #     return self

    @staticmethod
    def score_orthology_group(orthology_group, orthology_table):
        return orthology_group.assign_score(orthology_table)

    def rank_groups_by_score(self):  # TODO parallelize
        with ThreadPoolExecutor() as executor:
            executor.map(OrthologyGroups.score_orthology_group, self.orthology_groups_list)
<<<<<<< HEAD
        self.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)  # TODO does not work score is None
        return self
=======
        self.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)
        return self.orthology_groups_list
>>>>>>> 2af5a9c7ae422791ae65da0edd6a58ee9a6a48d2

    def remove_full_groups(self):
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(group.is_full): group for group in self}

            for future in as_completed(futures):
                group = futures[future]
                if future.result():
                    self.orthology_groups_list.remove(group)

        return self

    @staticmethod
    def create_custom_blast_db(orthology_group):
        orthology_group.fasta_to_custom_blast_db()

    def create_custom_blast_dbs(self):
        with ThreadPoolExecutor() as executor:
            executor.map(OrthologyGroups.create_custom_blast_db, self.orthology_groups_list)

    @staticmethod
    def write_ortholog_sequences_to_fasta(orthology_group):
        orthology_group.write_ortholog_sequences_to_fasta()

    def write_ortholog_sequences_to_fastas(self):
        with ThreadPoolExecutor() as executor:
            executor.map(OrthologyGroups.write_ortholog_sequences_to_fasta, self.orthology_groups_list)

