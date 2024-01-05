import pandas as pd


class GSE11223_MetadataProcessor:
    def process(self, metadata: pd.DataFrame) -> pd.DataFrame:
        SAMPLE_ID = 'geo_accession'
        PATIENT_ID = 'characteristics_ch1.0.patient'
        DISEASE = 'characteristics_ch1.51.disease'
        RESPONSE = None
        TIME_OF_BIOPSY = None

        metadata = metadata[[SAMPLE_ID, PATIENT_ID, DISEASE]]
        metadata['response'] = None
        metadata['time_of_biopsy'] = None
        
        metadata = metadata.rename(columns={
            SAMPLE_ID: 'sample_id',
            PATIENT_ID: 'patient_id',
            DISEASE: 'disease'
        })

        return metadata
