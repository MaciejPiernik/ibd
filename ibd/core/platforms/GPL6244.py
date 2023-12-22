
class GPL6244():
    def get_entrez_id_lists(self, locations):
        ids = locations.gene_assignment.map(self._gene_assignment_row_to_ids)

        return ids

    def _gene_assignment_row_to_ids(self, row):
        gene_assignments = str(row).split(' /// ')

        ids = []
        for gene_assignment in gene_assignments:
            gene_assignment_entries = gene_assignment.split(' // ')
            if len(gene_assignment_entries) == 5:
                if gene_assignment_entries[4].isdigit():
                    ids.append(gene_assignment_entries[4])
        
        return list(set(ids))