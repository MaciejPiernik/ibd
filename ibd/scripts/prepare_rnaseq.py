"""
Build the RNA-seq validation cohorts.

Downloads the nine datasets selected by screen_rnaseq.py and writes each one in
the same form as the microarray datasets in db/: a samples x genes parquet on a
log2 scale plus a metadata CSV giving disease and inflammation status. The
meta-analysis pipeline then runs on db_rnaseq/ without any modification.

Every dataset deposits its data differently, so there is one function per
dataset, each readable on its own. They share only the four helpers below.

Run:
    python ibd/scripts/prepare_rnaseq.py
"""

import os
import re
import tarfile
import urllib.request

import numpy as np
import pandas as pd

RAW = 'db_rnaseq_raw'
OUT = 'db_rnaseq'
HEADERS = 'data/geo_rnaseq_headers'
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series'


def download(gse, filename):
    """Fetch a supplementary file once, return the local path."""
    path = os.path.join(RAW, filename)
    if not os.path.exists(path):
        prefix = re.sub(r'\d{1,3}$', 'nnn', gse)
        print(f'  downloading {filename}')
        urllib.request.urlretrieve(f'{GEO}/{prefix}/{gse}/suppl/{filename}', path)

    return path


def series_field(gse, field):
    """One row of the series matrix header, as a list with one entry per sample."""
    for line in open(f'{HEADERS}/{gse}.txt', encoding='utf-8'):
        if line.startswith(f'!Sample_{field}'):
            return [value.strip('"') for value in line.rstrip('\n').split('\t')[1:]]


def series_characteristic(gse, name):
    """One characteristics field of the series matrix, one entry per sample."""
    for line in open(f'{HEADERS}/{gse}.txt', encoding='utf-8'):
        if line.startswith('!Sample_characteristics_ch1') and f'"{name}:' in line:
            values = [value.strip('"') for value in line.rstrip('\n').split('\t')[1:]]
            return [value.split(':', 1)[1].strip() for value in values]


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


def normalise(matrix, counts=True):
    """genes x samples -> samples x genes, log2 scale, low-expression genes dropped."""
    if counts:
        matrix = matrix / matrix.sum() * 1e6
    matrix = matrix.groupby(level=0).mean()
    matrix = matrix[(matrix >= 1).mean(axis=1) >= 0.2]

    return np.log2(matrix + 1).T


def save(gse, expression, disease, inflammation):
    expression.index = [f'{gse}_{sample}' for sample in expression.index]
    pd.DataFrame({
        'patient_id': expression.index,
        'disease': disease,
        'treatment': None,
        'response': None,
        'time_of_biopsy': None,
        'tissue': 'Colon',
        'inflammation': inflammation,
    }, index=expression.index).to_csv(f'{OUT}/{gse}_metadata.csv')
    expression.to_parquet(f'{OUT}/{gse}.parquet')

    counts = pd.Series(list(zip(disease, inflammation))).value_counts().to_dict()
    print(f'  {gse}: {expression.shape[0]} samples, {expression.shape[1]} genes, {counts}')


def per_sample_tpm(gse, tar_name):
    """Datasets depositing one TPM file per sample inside a RAW tar."""
    columns = {}
    with tarfile.open(download(gse, tar_name)) as tar:
        for member in tar.getmembers():
            values = pd.read_csv(tar.extractfile(member), sep='\t',
                                 index_col=0, compression='gzip').iloc[:, 0]
            columns[member.name.split('_')[0]] = values

    return pd.DataFrame(columns)


# ── the nine cohorts ──────────────────────────────────────────────────────

def gse109142():
    """PROTECT: treatment-naive paediatric UC, rectal biopsies. All UC active."""
    tpm = per_sample_tpm('GSE109142', 'GSE109142_RAW.tar')
    gsms = series_field('GSE109142', 'geo_accession')
    is_control = ['Ctl' in s for s in series_field('GSE109142', 'source_name_ch1')]
    control = dict(zip(gsms, is_control))

    expression = normalise(tpm, counts=False)
    disease = ['Ctrl' if control[s] else 'UC' for s in expression.index]
    save('GSE109142', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse117993():
    """Treatment-naive paediatric IBD, rectal biopsies at diagnostic colonoscopy."""
    tpm = per_sample_tpm('GSE117993', 'GSE117993_RAW.tar')
    diagnosis = dict(zip(series_field('GSE117993', 'geo_accession'),
                         series_characteristic('GSE117993', 'diagnosis')))

    keep = [s for s in tpm.columns if diagnosis[s] in ('UC', 'Control')]
    expression = normalise(tpm[keep], counts=False)
    disease = ['UC' if diagnosis[s] == 'UC' else 'Ctrl' for s in expression.index]
    save('GSE117993', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse165512():
    """Colonic biopsies from UC, CD and healthy donors."""
    matrix = pd.read_csv(download('GSE165512', 'GSE165512_normalized_counts.tsv.gz'),
                         sep='\t', index_col=0).drop(columns='Geneid')
    source = dict(zip(series_field('GSE165512', 'title'),
                      series_field('GSE165512', 'source_name_ch1')))

    keep = [s for s in matrix.columns if source[s] in ('colon_UC', 'colon_healthy')]
    expression = normalise(matrix[keep], counts=False)
    disease = ['UC' if source[s] == 'colon_UC' else 'Ctrl' for s in expression.index]
    save('GSE165512', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse193677():
    """MSCCR: adult registry, colonic biopsies, inflamed and non-inflamed UC."""
    counts = pd.read_csv(download('GSE193677', 'GSE193677_MSCCR_Biopsy_counts.txt.gz'),
                         sep=r'\s+', index_col=0)
    annotation = pd.DataFrame({
        'sample': [t.split(',')[0].strip() for t in series_field('GSE193677', 'title')],
        'group': series_characteristic('GSE193677', 'diseasetypere'),
        'region': series_characteristic('GSE193677', 'regionre'),
    }).set_index('sample')

    groups = {'UC.I': ('UC', 'Inflamed'), 'UC.NonI': ('UC', 'Uninflamed'),
              'Control.NonI': ('Ctrl', 'Uninflamed')}
    annotation = annotation[annotation['group'].isin(groups) & (annotation['region'] != 'Ileum')]

    # One biopsy per participant and group; several were sampled repeatedly.
    order = {'Rectum': 0, 'Sigmoid': 1, 'LeftColon': 2, 'Transverse': 3, 'RightColon': 4, 'Cecum': 5}
    annotation['participant'] = annotation.index.str.extract(r'^(MSCCR_reGRID_\d+)', expand=False)
    annotation = (annotation.assign(rank=annotation['region'].map(order))
                  .sort_values('rank').groupby(['participant', 'group']).head(1))

    expression = normalise(to_symbols(counts[list(annotation.index)]))
    labels = annotation.loc[expression.index, 'group']
    save('GSE193677', expression,
         [groups[g][0] for g in labels], [groups[g][1] for g in labels])


def gse128682():
    """Active UC, UC in remission and controls; colonic biopsies."""
    counts = pd.read_csv(download('GSE128682', 'GSE128682_raw_counts.tab.gz'),
                         sep='\t', index_col=0)
    expression = normalise(counts)

    groups = {'N': ('Ctrl', 'Uninflamed'), 'UC': ('UC', 'Inflamed'), 'R': ('UC', 'Uninflamed')}
    prefixes = [s.split('_')[0] for s in expression.index]
    save('GSE128682', expression,
         [groups[p][0] for p in prefixes], [groups[p][1] for p in prefixes])


def gse224758():
    """Active UC and healthy volunteers; colonic biopsies.

    The count matrix columns carry internal identifiers, so they are matched to
    the sample titles by position; the identifiers themselves corroborate it.
    """
    counts = pd.read_csv(download('GSE224758', 'GSE224758_read_counts.csv.gz'), index_col=0)
    titles = series_field('GSE224758', 'title')
    is_uc = dict(zip(counts.columns, ['ulcerative' in t for t in titles]))

    expression = normalise(to_symbols(counts))
    disease = ['UC' if is_uc[s] else 'Ctrl' for s in expression.index]
    save('GSE224758', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse334803():
    """UC and healthy controls; colonic biopsies.

    The deposited matrix carries two extra control libraries that are not part of
    the series, so the columns are restricted to the declared library names.
    """
    matrix = pd.read_excel(download('GSE334803', 'GSE334803_UCHC.xlsx'))
    counts = matrix.set_index('gene_name').filter(like='count.')
    counts.columns = [c.replace('count.', '') for c in counts.columns]
    libraries = [d.replace('Library name:', '').strip()
                 for d in series_field('GSE334803', 'description')]

    expression = normalise(counts[libraries])
    disease = ['UC' if s.startswith('UC') else 'Ctrl' for s in expression.index]
    save('GSE334803', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse266325():
    """Treatment-naive paediatric IBD, rectal biopsies. Titles give 'column, diagnosis'."""
    matrix = pd.read_csv(download('GSE266325', 'GSE266325_raw_counts.txt.gz'), sep='\t')
    diagnosis = dict(t.split(',') for t in series_field('GSE266325', 'title'))
    diagnosis = {k.strip(): v.strip() for k, v in diagnosis.items()}

    counts = matrix.set_index('gene_name')[[c for c in diagnosis if c in matrix.columns]]
    keep = [c for c in counts.columns if diagnosis[c] in ('UC', 'non IBD')]

    expression = normalise(counts[keep])
    disease = ['UC' if diagnosis[s] == 'UC' else 'Ctrl' for s in expression.index]
    save('GSE266325', expression, disease,
         ['Uninflamed' if d == 'Ctrl' else 'Inflamed' for d in disease])


def gse164918():
    """Quiescent UC and controls, sigmoid colon.

    A wound-healing study; only the 0 h biopsies precede the intervention, and
    the UC patients were in remission, so they form an uninflamed UC arm.
    """
    matrix = pd.read_csv(download('GSE164918', 'GSE164918_feaureCounts_RNAseq_matrix.txt.gz'),
                         sep='\t', comment='#', index_col=0)
    counts = matrix.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
    counts = counts[[c for c in counts.columns if '_0h_' in c]]

    expression = normalise(counts)
    disease = ['Ctrl' if s.startswith('Control') else 'UC' for s in expression.index]
    save('GSE164918', expression, disease, ['Uninflamed'] * len(disease))


if __name__ == '__main__':
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    for build in [gse109142, gse117993, gse165512, gse193677, gse128682,
                  gse224758, gse334803, gse266325, gse164918]:
        build()
