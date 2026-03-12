import mygene


class GPL5760():
    def get_entrez_id_lists(self, locations):
        mg = mygene.MyGeneInfo()

        accessions = locations['GB_ACC'].tolist()
        results = mg.querymany(accessions, scopes='refseq', fields='entrezgene', species='human')

        acc_to_entrez = {r['query']: [r.get('entrezgene', None)] for r in results}
        
        ids = locations['GB_ACC'].map(acc_to_entrez)

        return ids