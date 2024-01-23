import pandas as pd


class GPL10558():
    def get_entrez_id_lists(self, locations):
        ids = locations.Entrez_Gene_ID.map(lambda x: [str(int(x))] if not pd.isnull(x) else [None])

        return ids