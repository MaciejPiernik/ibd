"""
Run the meta-analysis pipeline over the RNA-seq validation cohorts.

The pipeline is used exactly as it was for the microarray analyses. This script
only selects the two comparisons and sets the minimum-dataset threshold to 50%
of the contributing cohorts, the same rule applied in the discovery analyses.

Run:
    RPY2_CFFI_MODE=ABI python ibd/scripts/run_validation.py
"""

import experiments.Experiment as exp
from experiments.Experiment import build_config

# output directory, inflammation status of the UC arm, number of cohorts
COMPARISONS = [
    ('rnaseq_uc_vs_hc', 'Inflamed', 8),
    ('rnaseq_uninf_vs_hc', 'Uninflamed', 3),
]

for output, inflammation, cohorts in COMPARISONS:
    exp.GROUP_1, exp.INFLAMMATION_1 = 'UC', inflammation
    exp.GROUP_2, exp.INFLAMMATION_2 = 'Ctrl', 'Uninflamed'

    exp.experiment(build_config(
        data_dir='./db_rnaseq',
        output_dir=f'results/revision_v2/{output}',
        skip_gene_mapping=True,
        min_datasets=-(-cohorts // 2),
    ))
