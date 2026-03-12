import pandas as pd


class GPL14951():
    def get_entrez_id_lists(self, locations):
        ids = locations.Entrez_Gene_ID.map(lambda x: [None] if pd.isna(x) else [int(x)])

        return ids