import tarfile

import numpy as np
import pandas as pd


MIN_EXPRESSION = 1          # counts per million, or transcripts per million
MIN_FRACTION_EXPRESSED = 0.2


def normalise(matrix, counts=True):
    """genes x samples -> samples x genes, log2 scale, low-expression genes dropped.

    Counts are converted to counts per million first; datasets deposited as
    transcripts per million are already on a comparable scale.
    """
    if counts:
        matrix = matrix / matrix.sum() * 1e6
    matrix = matrix.groupby(level=0).mean()
    matrix = matrix[(matrix >= MIN_EXPRESSION).mean(axis=1) >= MIN_FRACTION_EXPRESSED]

    return np.log2(matrix + 1).T


def to_symbols(matrix):
    """Replace Ensembl gene IDs in the index with gene symbols."""
    import mygene

    ids = [str(i).split('.')[0] for i in matrix.index]
    info = mygene.MyGeneInfo().querymany(
        ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
    symbols = info['symbol'].dropna()
    symbols = symbols[~symbols.index.duplicated()]
    matrix.index = pd.Index(ids).map(symbols)

    return matrix[matrix.index.notna()]


def build(dataset, expression, disease, inflammation):
    """Prefix the sample ids and assemble the harmonised metadata table."""
    expression.index = [f'{dataset.id}_{sample}' for sample in expression.index]

    metadata = pd.DataFrame({
        'patient_id': expression.index,
        'disease': disease,
        'treatment': None,
        'response': None,
        'time_of_biopsy': None,
        'tissue': 'Colon',
        'inflammation': inflammation,
    }, index=expression.index)

    return expression, metadata


def inflamed_disease(disease):
    """Controls are uninflamed by definition; disease samples are the inflamed arm."""
    return ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease]


def per_sample_tpm(dataset, tar_name):
    """Datasets depositing one transcripts-per-million file per sample in a tar."""
    columns = {}
    with tarfile.open(dataset.supplementary(tar_name)) as tar:
        for member in tar.getmembers():
            values = pd.read_csv(tar.extractfile(member), sep='\t',
                                 index_col=0, compression='gzip').iloc[:, 0]
            columns[member.name.split('_')[0]] = values

    return pd.DataFrame(columns)


class GSE109142_RnaSeqProcessor:
    """PROTECT: treatment-naive paediatric UC, rectal biopsies. All UC is active."""
    def process(self, dataset):
        tpm = per_sample_tpm(dataset, 'GSE109142_RAW.tar')
        control = dict(zip(dataset.sample_field('geo_accession'),
                           ['Ctl' in s for s in dataset.sample_field('source_name_ch1')]))

        expression = normalise(tpm, counts=False)
        disease = ['Ctrl' if control[s] else 'UC' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE117993_RnaSeqProcessor:
    """Treatment-naive paediatric IBD, rectal biopsies at diagnostic colonoscopy."""
    def process(self, dataset):
        tpm = per_sample_tpm(dataset, 'GSE117993_RAW.tar')
        diagnosis = dict(zip(dataset.sample_field('geo_accession'),
                             dataset.sample_characteristic('diagnosis')))

        keep = [s for s in tpm.columns if diagnosis[s] in ('UC', 'Control')]
        expression = normalise(tpm[keep], counts=False)
        disease = ['UC' if diagnosis[s] == 'UC' else 'Ctrl' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE165512_RnaSeqProcessor:
    """Colonic biopsies from UC, CD and healthy donors; ileal samples excluded."""
    def process(self, dataset):
        matrix = pd.read_csv(dataset.supplementary('GSE165512_normalized_counts.tsv.gz'),
                             sep='\t', index_col=0).drop(columns='Geneid')
        source = dict(zip(dataset.sample_field('title'), dataset.sample_field('source_name_ch1')))

        keep = [s for s in matrix.columns if source[s] in ('colon_UC', 'colon_healthy')]
        expression = normalise(matrix[keep], counts=False)
        disease = ['UC' if source[s] == 'colon_UC' else 'Ctrl' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE193677_RnaSeqProcessor:
    """MSCCR: adult registry, colonic biopsies, inflamed and non-inflamed UC.

    Participants were sampled at several sites, so one biopsy is kept per
    participant and group, preferring the most distal region available.
    """
    GROUPS = {'UC.I': ('UC', 'Inflamed'), 'UC.NonI': ('UC', 'Uninflamed'),
              'Control.NonI': ('Ctrl', 'Uninflamed')}
    REGIONS = {'Rectum': 0, 'Sigmoid': 1, 'LeftColon': 2,
               'Transverse': 3, 'RightColon': 4, 'Cecum': 5}

    def process(self, dataset):
        counts = pd.read_csv(dataset.supplementary('GSE193677_MSCCR_Biopsy_counts.txt.gz'),
                             sep=r'\s+', index_col=0)
        annotation = pd.DataFrame({
            'sample': [t.split(',')[0].strip() for t in dataset.sample_field('title')],
            'group': dataset.sample_characteristic('diseasetypere'),
            'region': dataset.sample_characteristic('regionre'),
        }).set_index('sample')

        annotation = annotation[annotation['group'].isin(self.GROUPS)
                                & (annotation['region'] != 'Ileum')]
        annotation['participant'] = annotation.index.str.extract(r'^(MSCCR_reGRID_\d+)', expand=False)
        annotation = (annotation.assign(rank=annotation['region'].map(self.REGIONS))
                      .sort_values('rank').groupby(['participant', 'group']).head(1))

        expression = normalise(to_symbols(counts[list(annotation.index)]))
        groups = annotation.loc[expression.index, 'group']

        return build(dataset, expression,
                     [self.GROUPS[g][0] for g in groups], [self.GROUPS[g][1] for g in groups])


class GSE128682_RnaSeqProcessor:
    """Active UC, UC in remission and controls; colonic biopsies."""
    GROUPS = {'N': ('Ctrl', 'Uninflamed'), 'UC': ('UC', 'Inflamed'), 'R': ('UC', 'Uninflamed')}

    def process(self, dataset):
        counts = pd.read_csv(dataset.supplementary('GSE128682_raw_counts.tab.gz'),
                             sep='\t', index_col=0)
        expression = normalise(counts)
        prefixes = [s.split('_')[0] for s in expression.index]

        return build(dataset, expression,
                     [self.GROUPS[p][0] for p in prefixes], [self.GROUPS[p][1] for p in prefixes])


class GSE224758_RnaSeqProcessor:
    """Active UC and healthy volunteers; colonic biopsies.

    The count matrix columns carry internal identifiers, so they are matched to
    the sample titles by position; the identifiers themselves corroborate it.
    """
    def process(self, dataset):
        counts = pd.read_csv(dataset.supplementary('GSE224758_read_counts.csv.gz'), index_col=0)
        is_uc = dict(zip(counts.columns, ['ulcerative' in t for t in dataset.sample_field('title')]))

        expression = normalise(to_symbols(counts))
        disease = ['UC' if is_uc[s] else 'Ctrl' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE334803_RnaSeqProcessor:
    """UC and healthy controls; colonic biopsies.

    The deposited matrix carries two extra control libraries that are not part
    of the series, so the columns are restricted to the declared library names.
    """
    def process(self, dataset):
        matrix = pd.read_excel(dataset.supplementary('GSE334803_UCHC.xlsx'))
        counts = matrix.set_index('gene_name').filter(like='count.')
        counts.columns = [c.replace('count.', '') for c in counts.columns]
        libraries = [d.replace('Library name:', '').strip()
                     for d in dataset.sample_field('description')]

        expression = normalise(counts[libraries])
        disease = ['UC' if s.startswith('UC') else 'Ctrl' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE266325_RnaSeqProcessor:
    """Treatment-naive paediatric IBD, rectal biopsies.

    Sample titles are 'matrix column, diagnosis'; Crohn's disease is excluded.
    """
    def process(self, dataset):
        matrix = pd.read_csv(dataset.supplementary('GSE266325_raw_counts.txt.gz'), sep='\t')
        diagnosis = {column.strip(): value.strip()
                     for column, value in (t.split(',') for t in dataset.sample_field('title'))}

        keep = [c for c in diagnosis if c in matrix.columns and diagnosis[c] in ('UC', 'non IBD')]
        expression = normalise(matrix.set_index('gene_name')[keep])
        disease = ['UC' if diagnosis[s] == 'UC' else 'Ctrl' for s in expression.index]

        return build(dataset, expression, disease, inflamed_disease(disease))


class GSE164918_RnaSeqProcessor:
    """Quiescent UC and controls, sigmoid colon.

    A wound-healing study: only the 0 h biopsies precede the intervention, and
    the UC patients were in remission, so they form an uninflamed UC arm.
    """
    def process(self, dataset):
        matrix = pd.read_csv(dataset.supplementary('GSE164918_feaureCounts_RNAseq_matrix.txt.gz'),
                             sep='\t', comment='#', index_col=0)
        counts = matrix.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])

        expression = normalise(counts[[c for c in counts.columns if '_0h_' in c]])
        disease = ['Ctrl' if s.startswith('Control') else 'UC' for s in expression.index]

        return build(dataset, expression, disease, ['Uninflamed'] * len(disease))
