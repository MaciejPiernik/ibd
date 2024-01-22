import pandas as pd


class GPL13158():
    def get_entrez_id_lists(self, locations):
        ids = locations.ENTREZ_GENE_ID.map(lambda x: str(x).split(' /// '))

        return ids