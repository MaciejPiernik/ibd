"""
Check the numbers quoted in the manuscript against the analysis results.

Gene-level effect sizes are looked up directly and reported as matching or not.
Everything else numeric (enrichment scores, correlations, counts) is listed with
its line number and surrounding text, so the remaining claims can be worked
through against the results rather than by eye.

Run:
    python ibd/scripts/_check_manuscript_numbers.py --manuscript ../uc_paper/main.tex
"""

import argparse
import re

import pandas as pd

# Effect-size notation used in the text, and the analysis it refers to.
SOURCES = {
    'g': 'uc_vs_hc',
    'g_{\\mathrm{infl}}': 'uc_vs_hc',
    'g_{\\mathrm{uninfl}}': 'uninf_vs_hc',
}

EFFECT = re.compile(r'\$(g(?:_\{\\mathrm\{\w+\}\})?)\s*=\s*([+-]?\d+\.\d+)\$')
SYMBOL = re.compile(r'\b([A-Z][A-Z0-9]{2,9})\b')
OTHER = re.compile(r'NES|Pearson|Spearman|I\^2|q = |\d{1,3}(?:\{,\})?\d{3}|\b\d+\s*(?:datasets|pathways|genes|samples)\b')


def check_effects(text, results, tolerance):
    tables = {name: pd.read_csv(f'{results}/{name}/gene_meta_all.csv').set_index('gene')
              for name in set(SOURCES.values())}

    rows = []
    for match in EFFECT.finditer(text):
        analysis = SOURCES.get(match.group(1))
        if analysis is None:
            continue

        # the gene is the last symbol-like token before the value
        preceding = SYMBOL.findall(text[max(0, match.start() - 200):match.start()])
        gene = next((s for s in reversed(preceding) if s in tables[analysis].index), None)
        quoted = float(match.group(2))

        if gene is None:
            rows.append({'gene': '?', 'analysis': analysis, 'quoted': quoted,
                         'current': None, 'status': 'gene not identified'})
            continue

        current = tables[analysis].loc[gene, 'hedges_g']
        rows.append({'gene': gene, 'analysis': analysis, 'quoted': quoted,
                     'current': round(current, 3),
                     'status': 'ok' if abs(current - quoted) < tolerance
                               else f'differs by {abs(current - quoted):.2f}'})

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manuscript', default='../uc_paper/main.tex')
    parser.add_argument('--results', default='results/revision_v2')
    parser.add_argument('--tolerance', type=float, default=0.005)
    parser.add_argument('--inventory', action='store_true',
                        help='list every remaining numeric claim with its line number')
    args = parser.parse_args()

    text = open(args.manuscript, encoding='utf-8').read()

    effects = check_effects(text, args.results, args.tolerance)
    changed = effects[effects['status'] != 'ok']
    print(f'Gene effect sizes quoted: {len(effects)}  '
          f'(unchanged {int((effects["status"] == "ok").sum())}, '
          f'need updating {len(changed)})\n')
    if len(changed):
        print(changed.to_string(index=False))

    if args.inventory:
        print('\n\nRemaining numeric claims, by line:\n')
        for number, line in enumerate(text.split('\n'), 1):
            stripped = line.strip()
            if not stripped or stripped.startswith('%') or stripped.startswith('\\'):
                continue
            if OTHER.search(stripped) and not EFFECT.search(stripped):
                print(f'{number:5d}  {stripped[:150]}')


if __name__ == '__main__':
    main()
