import numpy as np
import pandas as pd


DISEASE_MAP = {
    'CD': 'CD',
    "Crohn's disease": "CD",

    'UC': 'UC',
    'ulcerative colitis': 'UC',
    'Ulcerative colitis': 'UC',
    'Ulcerative Colitis (UC)': 'UC',

    'Normal': 'Ctrl',
    'Control': 'Ctrl',
    'control': 'Ctrl',
    'Healthy control': 'Ctrl',
}


class GSE11223_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.0.patient'
        DISEASE = 'characteristics_ch1.51.disease'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID]
        result['disease'] = metadata[DISEASE].map(DISEASE_MAP)
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None

        return result


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

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID]
        result['disease'] = metadata[DISEASE].map(lambda x: 'UC' if 'UC' in x else 'Ctrl' if 'Control' in x else None)
        result['treatment'] = metadata[TREATMENT].map(lambda x: None if x == 'CO' else x)
        result['response'] = metadata[RESPONSE].map(lambda x: 'No' if 'NR' in x.replace(' ', '_').split('_') else ('Yes' if 'R' in x.replace(' ', '_').split('_') else ('Other' if 'other' in x.replace(' ', '_').split('_') else None)))
        result['time_of_biopsy'] = metadata[TIME_OF_BIOPSY].map(lambda x: None if x == 'CO' else x)

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


class GSE16879_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = 'characteristics_ch1.1.disease'
        TREATMENT = 'title'
        RESPONSE = 'characteristics_ch1.2.response to infliximab'
        TIME_OF_BIOPSY = 'characteristics_ch1.3.before or after first infliximab treatment'

        TIME_OF_BIOPSY_MAP = {
            'Before first infliximab treatment': 'Before',
            'After first infliximab treatment': 'After',
            'Not applicable': None
        }

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x.split('_')[0])
        result['disease'] = metadata[DISEASE].map(DISEASE_MAP)
        result['treatment'] = metadata[TREATMENT].map(lambda x: 'IFX' if len(x.split('_')) == 2 else None)
        result['response'] = metadata[RESPONSE].map(lambda x: 'Yes' if x == 'Yes' else ('No' if x == 'No' else None))
        result['time_of_biopsy'] = metadata[TIME_OF_BIOPSY].map(TIME_OF_BIOPSY_MAP)

        return result


class GSE52746_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.2.patient'
        DISEASE = 'title'
        TREATMENT = 'title'
        RESPONSE = 'title'
        TIME_OF_BIOPSY = 'title'

        result = pd.DataFrame()

        result['patient_id'] = metadata.apply(lambda x: x['title'].split('(')[1][:-1] if 'control' in x['title'] else x[PATIENT_ID], axis=1)
        result['disease'] = metadata[DISEASE].map(lambda x: 'CD' if 'CD' in x else 'Ctrl')
        result['treatment'] = metadata[TREATMENT].map(lambda x: None if 'control' in x else 'anti-TNF')
        result['response'] = metadata[RESPONSE].map(lambda x: 'Yes' if (' with ' in x and 'Inactive' in x) else ('No' if (' with ' in x and 'Active' in x) else None))
        result['time_of_biopsy'] = metadata[TIME_OF_BIOPSY].map(lambda x: 'Before' if 'without' in x else ('After' if ' with ' in x else None))

        return result


class GSE36807_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = 'characteristics_ch1.0.diagnosis'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        DISEASE_MAP = {
            'Healthy control': 'Ctrl',
            "Crohn's Disease": 'CD',
            'Ulcerative colitis': 'UC'
        }

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x.split('_')[0])
        result['disease'] = metadata[DISEASE].map(DISEASE_MAP)
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None

        return result


class GSE22619_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = 'title'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        PATIENT_ID_MAP = {
            'Healthy individual of twinpair': '1',
            'Diseased individual (ulcerative colitis) of twinpair': '2',
        }

        DISEASE_MAP = {
            'Healthy individual of twinpair': 'Ctrl',
            'Diseased individual (ulcerative colitis) of twinpair': 'UC',
        }

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: f"{x.split(' #')[1]}_{PATIENT_ID_MAP[x.split(' #')[0]]}")
        result['disease'] = metadata[DISEASE].map(lambda x: DISEASE_MAP[x.split(' #')[0]])
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None

        return result
    

class GSE9452_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = [c for c in metadata.columns if c.startswith('characteristics')]
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x.split('_')[-1])
        result['disease'] = metadata.apply(lambda x: 'UC' if 'Ulcerative colitis' in x[DISEASE][x[DISEASE].isna() == False].index[0] else 'Ctrl', axis=1)
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None

        return result


class GSE72780_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = None
        DISEASE = 'source_name_ch1'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None

        result = pd.DataFrame()

        result['disease'] = metadata[DISEASE].map(lambda x: 'CD' if '(CD)' in x else 'Ctrl')
        result['patient_id'] = range(len(result))
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None

        return result
    

class GSE179285_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = 'characteristics_ch1.1.diagnosis'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None
        TISSUE = 'characteristics_ch1.3.tissue'
        INFLAMMATION = 'characteristics_ch1.2.inflammation'

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x.split('_')[0])
        result['disease'] = metadata[DISEASE].map(DISEASE_MAP)
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None
        result['tissue'] = metadata[TISSUE]
        result['inflammation'] = metadata[INFLAMMATION]

        return result


class GSE75214_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'title'
        DISEASE = 'characteristics_ch1.1.disease'
        TREATMENT = None
        RESPONSE = None
        TIME_OF_BIOPSY = None
        TISSUE = 'characteristics_ch1.0.tissue'
        INFLAMMATION = 'characteristics_ch1.2.disease activity'

        INFLAMMATION_MAP = {
            'active': 'Inflamed',
            'inactive': 'Uninflamed',
            'normal': 'Uninflamed'
        }

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x.split('_')[-1])
        result['disease'] = metadata[DISEASE].map(DISEASE_MAP)
        result['treatment'] = None
        result['response'] = None
        result['time_of_biopsy'] = None
        result['tissue'] = metadata[TISSUE].map(lambda x: x.capitalize())
        result['inflammation'] = metadata[INFLAMMATION].map(INFLAMMATION_MAP)

        return result


class GSE92415_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        PATIENT_ID = 'characteristics_ch1.1.subject'
        DISEASE = 'characteristics_ch1.2.disease'
        TREATMENT = 'characteristics_ch1.4.treatment'
        RESPONSE = 'characteristics_ch1.6.wk6response'
        TIME_OF_BIOPSY = 'characteristics_ch1.5.visit'
        TISSUE = 'characteristics_ch1.0.tissue'
        INFLAMMATION = None
        MAYO_SCORE = 'characteristics_ch1.7.mayo score'

        curr_uid = 0
        def get_next_id():
            nonlocal curr_uid
            curr_uid += 1
            return f'C{curr_uid}'

        result = pd.DataFrame()

        result['patient_id'] = metadata[PATIENT_ID].map(lambda x: x if not pd.isna(x) else get_next_id())
        result['disease'] = metadata[DISEASE].map(lambda x: 'UC' if not pd.isna(x) else 'Ctrl')
        result['treatment'] = metadata[TREATMENT].map(lambda x: x.capitalize() if not pd.isna(x) else None)
        result['response'] = metadata[RESPONSE].map(lambda x: x if not pd.isna(x) else None)
        result['time_of_biopsy'] = metadata[TIME_OF_BIOPSY].map(lambda x: 'W0' if x == 'Week 0' else ('W6' if x == 'Week 6' else None))
        result['tissue'] = metadata[TISSUE].map(lambda x: 'Colon')
        result['inflammation'] = None
        result['mayo_score'] = metadata[MAYO_SCORE]

        return result