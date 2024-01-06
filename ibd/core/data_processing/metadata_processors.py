import pandas as pd


DISEASE_MAP = {
    'CD': 'CD',
    "Crohn's disease": "CD",

    'UC': 'UC',
    'ulcerative colitis': "UC",

    'control': "Ctrl",
}


class GSE11223_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.0.patient'
        DISEASE = 'characteristics_ch1.51.disease'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        metadata = metadata[[PATIENT_ID, DISEASE]]
        metadata['treatment'] = None
        metadata['response'] = None
        metadata['time_of_biopsy'] = None
        
        metadata = metadata.rename(columns={
            PATIENT_ID: 'patient_id',
            DISEASE: 'disease'
        })

        return metadata

class GSE75214_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = None
        DISEASE = 'characteristics_ch1.1.disease'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        metadata = metadata[[DISEASE]]
        metadata['patient_id'] = range(len(metadata))
        metadata['treatment'] = None
        metadata['response'] = None
        metadata['time_of_biopsy'] = None
        
        metadata = metadata.rename(columns={
            DISEASE: 'disease'
        })

        metadata['disease'] = metadata['disease'].map(DISEASE_MAP)

        return metadata
    
class GSE73661_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.0.study individual number'
        DISEASE = 'title'
        TREATMENT = 'characteristics_ch1.2.induction therapy_maintenance therapy'
        RESPONSE = 'title'
        TIME_OF_BIOPSY = 'characteristics_ch1.1.week (w)'

        result = metadata[[PATIENT_ID, TREATMENT, TIME_OF_BIOPSY]]

        result = result.rename(columns={
            PATIENT_ID: 'patient_id',
            TREATMENT: 'treatment',
            TIME_OF_BIOPSY: 'time_of_biopsy'
        })

        result['disease'] = metadata[DISEASE].map(lambda x: 'UC' if 'UC' in x else 'Ctrl' if 'Control' in x else None)
        result['response'] = metadata[RESPONSE].map(lambda x: 'NR' if 'NR' in x.replace(' ', '_').split('_') else ('R' if 'R' in x.replace(' ', '_').split('_') else ('Other' if 'other' in x.replace(' ', '_').split('_') else None)))

        return result

class GSE23597_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.3.subject'
        DISEASE = None
        TREATMENT = 'characteristics_ch1.1.dose'
        RESPONSE = ['characteristics_ch1.4.wk8 response', 'characteristics_ch1.5.wk30 response']
        TIME_OF_BIOPSY = 'characteristics_ch1.2.time'

        result = metadata[[PATIENT_ID, TREATMENT, TIME_OF_BIOPSY]]

        result = result.rename(columns={
            PATIENT_ID: 'patient_id',
            TREATMENT: 'treatment',
            TIME_OF_BIOPSY: 'time_of_biopsy'
        })

        result['disease'] = 'UC'
        result['treatment'] = metadata[TREATMENT].map(lambda x: 'Placebo' if x == 'placebo' else 'IFX')
        result['response'] = metadata.apply(lambda x: x[RESPONSE[0]] if x[TIME_OF_BIOPSY] == 'W8' else (x[RESPONSE[1]] if x[TIME_OF_BIOPSY] == 'W30' else None), axis=1)

        return result