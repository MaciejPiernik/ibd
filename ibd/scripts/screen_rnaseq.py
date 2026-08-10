"""
Screen GEO for the RNA-seq validation cohorts.

The discovery analysis enumerated the microarray corpus exhaustively. The same
search and the same inclusion criteria are applied here to RNA-seq, with only
the platform restriction lifted, so the validation set is every qualifying
dataset rather than a selection from among them.

Criteria, inherited from the discovery analysis:
  - colonic or rectal mucosal biopsies
  - gene-level quantification deposited (discovery: probe-to-gene annotation)
  - samples unambiguously classifiable by disease status from the metadata
  - not a subset, combination or duplicate of another included study
  - bulk tissue, since the discovery analysis is bulk
  - at least 3 samples per group and 6 in total

Stage 1 excludes on the study description. Stage 2 required reading the sample
metadata and the deposited files, and those decisions are recorded explicitly in
DECISIONS below so the screen can be audited.

Writes the screening table and caches each candidate's series matrix header,
which prepare_rnaseq.py then reads.

Run:
    python ibd/scripts/screen_rnaseq.py
"""

import gzip
import io
import json
import os
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

QUERY = ('("ulcerative colitis"[All Fields] OR "inflammatory bowel disease"[All Fields]) '
         'AND "Homo sapiens"[Organism] '
         'AND "expression profiling by high throughput sequencing"[DataSet Type] '
         'AND gse[Entry Type]')
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series'
HEADERS = 'data/geo_rnaseq_headers'
OUTPUT = 'results/revision_v2'

# First matching rule is the recorded reason for exclusion.
STAGE1 = [
    ('single-cell or single-nucleus', r'single[ -]cell|single[ -]nucleus|scrna|snrna|10x genomics|cite-seq|smart-seq'),
    ('spatial transcriptomics', r'spatial transcriptom|visium|geomx'),
    ('not primary tissue', r'organoid|enteroid|colonoid|cell line|caco-2|ht-29|hct116|in vitro|cultured cell'),
    ('sorted or isolated cells', r'sorted|facs|flow[- ]sorted|isolated (immune|epithelial|t )|cd4\+|cd8\+'),
    ('non-coding RNA profiling', r'\bmirna\b|micrornas?\b|lncrna profiling|small rna'),
]

# Datasets reaching stage 2, with the decision and the criterion behind it.
DECISIONS = {
    'GSE109142': ('include', 'Treatment-naive paediatric UC, rectal biopsies'),
    'GSE117993': ('include', 'Treatment-naive paediatric IBD, rectal biopsies'),
    'GSE165512': ('include', 'Colonic biopsies, UC and healthy'),
    'GSE193677': ('include', 'Adult registry, colonic biopsies, inflamed and non-inflamed UC'),
    'GSE128682': ('include', 'Colonic biopsies, active UC, UC in remission and controls'),
    'GSE224758': ('include', 'Colonic biopsies, active UC and healthy volunteers'),
    'GSE334803': ('include', 'Colonic biopsies, UC and healthy controls'),
    'GSE266325': ('include', 'Treatment-naive paediatric IBD, rectal biopsies'),
    'GSE164918': ('include', 'Quiescent UC and controls, sigmoid colon; baseline biopsies only'),
    'GSE95437': ('exclude', 'No gene-level quantification deposited (CAGE tag coordinates only)'),
    'GSE296485': ('exclude', 'No expression matrix deposited (differential expression list only)'),
    'GSE83687': ('exclude', 'Surgical resection specimens, not mucosal biopsies'),
    'GSE166925': ('exclude', 'Surgical resection specimens, not mucosal biopsies'),
    'GSE57945': ('exclude', 'Ileal biopsies only'),
    'GSE137344': ('exclude', 'Ileal biopsies only'),
    'GSE81266': ('exclude', 'Ileal pouch tissue, not colonic mucosa'),
    'GSE135223': ('exclude', 'Isolated macrophages, not whole mucosa'),
    'GSE277415': ('exclude', 'Sorted epithelial and lamina propria populations'),
}


def search():
    """Every GEO series matching the discovery search, restricted to RNA-seq."""
    url = f'{EUTILS}/esearch.fcgi?db=gds&retmax=500&term={urllib.parse.quote(QUERY)}'
    ids = re.findall(r'<Id>(\d+)</Id>', urllib.request.urlopen(url, timeout=120).read().decode())

    records = []
    for start in range(0, len(ids), 100):
        batch = ','.join(ids[start:start + 100])
        summary = json.load(urllib.request.urlopen(
            f'{EUTILS}/esummary.fcgi?db=gds&retmode=json&id={batch}', timeout=120))
        records += [summary['result'][uid] for uid in summary['result']['uids']]
        time.sleep(0.5)

    return records


def series_header(gse):
    """Cache the !Sample_ header of a series matrix, without the expression table."""
    path = f'{HEADERS}/{gse}.txt'
    if os.path.exists(path):
        return

    prefix = re.sub(r'\d{1,3}$', 'nnn', gse)
    request = urllib.request.Request(
        f'{GEO}/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz',
        headers={'Accept-Encoding': 'identity', 'User-Agent': 'ibd-meta-analysis/1.0'})

    raw = io.BytesIO()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                while True:
                    chunk = response.read(65536)
                    if not chunk:
                        break
                    raw.write(chunk)
                    if b'!series_matrix_table_begin' in raw.getvalue()[-200000:]:
                        break
            break
        except urllib.error.HTTPError:      # NCBI throttles bursts of requests
            time.sleep(5 * (attempt + 1))
            raw = io.BytesIO()

    raw.seek(0)
    text = gzip.zlib.decompressobj(gzip.zlib.MAX_WBITS | 16).decompress(raw.read()).decode(
        'utf-8', errors='replace')
    open(path, 'w', encoding='utf-8').write(text.split('!series_matrix_table_begin')[0])


if __name__ == '__main__':
    os.makedirs(HEADERS, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    records = search()
    print(f'GEO returned {len(records)} series')

    rows = []
    for record in records:
        gse = record['accession']
        text = f"{record['title']} {record['summary']}".lower()
        reason = next((name for name, pattern in STAGE1 if re.search(pattern, text)), None)
        if not reason and record['n_samples'] < 6:
            reason = 'fewer than 6 samples'
        if not reason and not re.search(r'ulcerative colitis|\buc\b|\bibd\b|inflammatory bowel', text):
            reason = 'no ulcerative colitis samples'
        if not reason and not re.search(r'biops|mucosa|colon|rect|intestin|sigmoid|ileum|gut', text):
            reason = 'not intestinal mucosa'

        decision, note = DECISIONS.get(gse, ('exclude', reason or 'no UC vs control colonic comparison'))
        rows.append({'accession': gse, 'n_samples': record['n_samples'], 'decision': decision,
                     'reason': note if gse in DECISIONS else reason or note,
                     'title': record['title'][:120]})

    table = pd.DataFrame(rows).sort_values(['decision', 'n_samples'], ascending=[True, False])
    table.to_csv(f'{OUTPUT}/rnaseq_screening.csv', index=False)

    included = table[table['decision'] == 'include']
    print(f'Included: {len(included)}')
    print(included[['accession', 'n_samples', 'reason']].to_string(index=False))

    for gse in table['accession']:
        if gse in DECISIONS:
            series_header(gse)
    print(f'\nSeries headers cached in {HEADERS}/')
