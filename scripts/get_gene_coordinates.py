import gzip, csv

ids = {}
with open('data/reference/gene_symbol_to_id.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ids[row['ensembl_gene_id']] = row['gene_symbol']

results = []
with gzip.open('data/reference/Bos_taurus.ARS-UCD1.2.110.gtf.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] != 'gene':
            continue
        attrs = fields[8]
        for gid in ids:
            if f'gene_id "{gid}"' in attrs:
                results.append({
                    'gene_symbol': ids[gid],
                    'ensembl_gene_id': gid,
                    'chromosome': fields[0],
                    'start': fields[3],
                    'end': fields[4],
                    'strand': fields[6]
                })

with open('data/reference/candidate_gene_coordinates_provisional.tsv', 'w') as out:
    out.write("gene_symbol\tensembl_gene_id\tchromosome\tstart\tend\tstrand\n")
    for r in sorted(results, key=lambda x: x['gene_symbol']):
        out.write(f"{r['gene_symbol']}\t{r['ensembl_gene_id']}\t{r['chromosome']}\t{r['start']}\t{r['end']}\t{r['strand']}\n")

print(open('data/reference/candidate_gene_coordinates_provisional.tsv').read())
print(f"\nFound {len(results)} of 13 genes in the GTF")
