#!/usr/bin/env python3
"""H9: Candidate-level pairwise breed Weir-Cockerham FST, 16 breeds."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import math
import subprocess
import tempfile

import pandas as pd

from H7_candidate_regional_fst import ROOT, INPUT, VCFS, EXPECTED, samples, extract, stats, independent_wc
from H8_candidate_country_fst import load_manifest

OUT = ROOT / 'results/H_integration/H9_candidate_breed_fst'


def main():
    candidates = pd.read_csv(INPUT, sep='\t', dtype={'CHROM': str})
    required = {'candidate_id', 'SYMBOL', 'CHROM', 'POS', 'REF', 'ALT'}
    if len(candidates) != 62 or not required.issubset(candidates.columns):
        raise ValueError('Candidate table must contain 62 variants and required fields')
    if candidates.duplicated(['CHROM', 'POS', 'REF', 'ALT']).any():
        raise ValueError('Duplicate candidate variants')
    if candidates['candidate_id'].duplicated().any():
        raise ValueError('Duplicate candidate IDs')

    manifest = load_manifest()
    if manifest['breed'].isna().any() or manifest['cohort'].isna().any():
        raise ValueError('Missing breed or cohort')
    breed_counts = manifest.groupby('breed').size().sort_index()
    breeds = list(breed_counts.index)
    if len(breeds) != 16 or len(manifest) != 346:
        raise ValueError('Expected 16 breeds and 346 animals')
    breed_cohorts = manifest.groupby('breed')['cohort'].nunique()
    if (breed_cohorts != 1).any():
        raise ValueError('Breed spans multiple cohort VCF sets; needs explicit handling')
    cohort_for = manifest.groupby('breed')['cohort'].first().to_dict()

    cohort_ids, paths = {}, {}
    for cohort, expected_n in EXPECTED.items():
        for gene in sorted(candidates['SYMBOL'].unique()):
            path = VCFS / cohort / f'{gene}_qc.vcf.gz'
            if not path.is_file():
                raise FileNotFoundError(path)
            ids = samples(path)
            if len(ids) != expected_n or len(set(ids)) != expected_n:
                raise ValueError(f'Bad sample list: {path}')
            if cohort in cohort_ids and ids != cohort_ids[cohort]:
                raise ValueError(f'VCF sample order differs: {path}')
            cohort_ids.setdefault(cohort, ids)
            paths[(cohort, gene)] = path
    if set(cohort_ids['europe']) & set(cohort_ids['south_africa']):
        raise ValueError('Overlapping cohort sample IDs')
    if set().union(*(set(v) for v in cohort_ids.values())) != set(manifest['vcf_id']):
        raise ValueError('Manifest IDs do not match VCF sample IDs')

    breed_ids = {}
    for breed in breeds:
        assigned = set(manifest.loc[manifest['breed'] == breed, 'vcf_id'])
        breed_ids[breed] = [sid for sid in cohort_ids[cohort_for[breed]] if sid in assigned]
        if len(breed_ids[breed]) != int(breed_counts[breed]):
            raise ValueError(f'Breed membership mismatch: {breed}')

    extracted = {}
    for row in candidates.itertuples(index=False):
        cid = str(row.candidate_id)
        extracted[cid] = {}
        for cohort, n in EXPECTED.items():
            calls = extract(paths[(cohort, row.SYMBOL)], str(row.CHROM), int(row.POS),
                            str(row.REF), str(row.ALT), n)
            extracted[cid][cohort] = None if calls is None else dict(zip(cohort_ids[cohort], calls))

    rows = []
    with tempfile.TemporaryDirectory(prefix='H9_breed_fst_') as td:
        tmp = Path(td)
        for pair_i, (breed_a, breed_b) in enumerate(combinations(breeds, 2), 1):
            ids_a, ids_b = breed_ids[breed_a], breed_ids[breed_b]
            if set(ids_a) & set(ids_b):
                raise ValueError('Overlapping breed groups')
            pair = f'{breed_a}__{breed_b}'
            (tmp / 'a.samples').write_text('\n'.join(ids_a) + '\n')
            (tmp / 'b.samples').write_text('\n'.join(ids_b) + '\n')
            eligible = []
            for candidate in candidates.itertuples(index=False):
                cid = str(candidate.candidate_id)
                source_a = extracted[cid][cohort_for[breed_a]]
                source_b = extracted[cid][cohort_for[breed_b]]
                calls_a = None if source_a is None else [source_a[sid] for sid in ids_a]
                calls_b = None if source_b is None else [source_b[sid] for sid in ids_b]
                sa, sb = stats(calls_a), stats(calls_b)
                result = {
                    'candidate_id': cid, 'SYMBOL': candidate.SYMBOL,
                    'CHROM': str(candidate.CHROM), 'POS': int(candidate.POS),
                    'REF': candidate.REF, 'ALT': candidate.ALT,
                    'breed_a': breed_a, 'breed_b': breed_b, 'comparison': pair,
                    'breed_a_n': len(ids_a), 'breed_b_n': len(ids_b),
                    'small_breed_pair': min(len(ids_a), len(ids_b)) < 5,
                    'breed_a_called': sa['called'], 'breed_b_called': sb['called'],
                    'breed_a_alt_af': sa['alt_af'], 'breed_b_alt_af': sb['alt_af'],
                    'breed_a_hom_ref': sa['hom_ref'], 'breed_a_het': sa['het'],
                    'breed_a_hom_alt': sa['hom_alt'],
                    'breed_b_hom_ref': sb['hom_ref'], 'breed_b_het': sb['het'],
                    'breed_b_hom_alt': sb['hom_alt'],
                    'fst_wc': None, 'fst_wc_independent': None,
                }
                if source_a is None:
                    result['fst_status'] = 'absent_breed_a_cohort_vcf'
                elif source_b is None:
                    result['fst_status'] = 'absent_breed_b_cohort_vcf'
                elif sa['called'] < 2 or sb['called'] < 2:
                    result['fst_status'] = 'insufficient_called_genotypes'
                else:
                    result['fst_status'] = 'pending'
                    eligible.append((result, calls_a + calls_b, sa, sb))
                rows.append(result)
            eligible.sort(key=lambda x: (int(x[0]['CHROM']) if x[0]['CHROM'].isdigit() else 10**9, x[0]['POS']))
            if eligible:
                vcf = tmp / 'pair.vcf'
                with vcf.open('w') as f:
                    f.write('##fileformat=VCFv4.2\n')
                    for chrom in sorted({x[0]['CHROM'] for x in eligible}):
                        f.write(f'##contig=<ID={chrom}>\n')
                    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Diploid genotype">\n')
                    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ids_a + ids_b) + '\n')
                    for r, calls, _, _ in eligible:
                        f.write('\t'.join([r['CHROM'], str(r['POS']), r['candidate_id'], r['REF'], r['ALT'], '.', 'PASS', '.', 'GT', *calls]) + '\n')
                prefix = tmp / 'pair_result'
                proc = subprocess.run(['vcftools', '--vcf', str(vcf),
                                       '--weir-fst-pop', str(tmp / 'a.samples'),
                                       '--weir-fst-pop', str(tmp / 'b.samples'),
                                       '--out', str(prefix)], capture_output=True, text=True)
                if proc.returncode:
                    raise RuntimeError(f'{pair}: {proc.stdout}\n{proc.stderr}')
                fst = pd.read_csv(str(prefix) + '.weir.fst', sep=r'\s+', dtype={'CHROM': str})
                lookup = {}
                for item in fst.itertuples(index=False):
                    key = (str(item.CHROM), int(item.POS))
                    if key in lookup:
                        raise ValueError(f'Duplicate VCFtools site: {pair} {key}')
                    lookup[key] = float(item.WEIR_AND_COCKERHAM_FST)
                expected_keys = {(r['CHROM'], r['POS']) for r, _, _, _ in eligible}
                if set(lookup) - expected_keys:
                    raise ValueError(f'Unexpected VCFtools positions: {pair}')
                for r, _, sa, sb in eligible:
                    key = (r['CHROM'], r['POS'])
                    value = lookup.get(key, float('nan'))
                    independent = independent_wc(sa, sb)
                    if not math.isfinite(value):
                        if math.isfinite(independent):
                            raise ValueError(f'VCFtools undefined but independent finite: {pair} {key}')
                        r['fst_status'] = 'not_estimable_vcftools'
                    else:
                        if not math.isfinite(independent) or abs(value - independent) > 1e-5:
                            raise ValueError(f'Estimator disagreement: {pair} {key}: {value} vs {independent}')
                        r['fst_wc'] = value
                        r['fst_wc_independent'] = independent
                        r['fst_status'] = 'estimated'
            if pair_i % 10 == 0 or pair_i == 120:
                print(f'Processed {pair_i}/120 breed pairs', flush=True)

    result = pd.DataFrame(rows)
    if len(result) != 62 * 120 or result.duplicated(['candidate_id', 'comparison']).any():
        raise ValueError('Incorrect row count or duplicate candidate-comparison')
    if (result['fst_status'] == 'pending').any():
        raise ValueError('Unresolved pending rows')
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / 'candidate_breed_pairwise_fst_62.tsv', sep='\t', index=False,
                  na_rep='NA', float_format='%.9f')
    manifest.to_csv(OUT / 'breed_sample_manifest_346.tsv', sep='\t', index=False)
    lines = [
        'H9: Candidate-level pairwise breed Weir-Cockerham FST',
        f'Samples: {len(manifest)}', f'Breeds: {len(breeds)}',
        'Candidate variants: 62', 'Breed pairs: 120',
        f'Candidate-pair records: {len(result)}', '', 'Status counts:',
        *[f'{k}: {v}' for k, v in result['fst_status'].value_counts().items()],
        '', f'Small-breed pair records (min breed n < 5): {int(result["small_breed_pair"].sum())}',
        'Every finite FST checked independently within 1e-5.',
        'Absent variant-only cohort VCF records were not imputed as reference.',
        'Negative estimates retained. Small-breed estimates require caution.',
        'Holstein-Friesian pooled as one breed across its three countries.',
    ]
    (OUT / 'candidate_breed_fst_run_summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n' + '\n'.join(lines))
    print('\nPASS: H9 completed')


if __name__ == '__main__':
    main()
