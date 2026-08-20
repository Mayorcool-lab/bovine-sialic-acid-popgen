import urllib.request, json, time

with open('data/reference/candidate_genes.txt') as f:
    genes = [line.strip() for line in f if line.strip()]

results = []
for gene in genes:
    url = f"https://rest.ensembl.org/xrefs/symbol/bos_taurus/{gene}?content-type=application/json"
    try:
        with urllib.request.urlopen(url) as resp:
            data = json.loads(resp.read())
        gene_ids = [d['id'] for d in data if d.get('type') == 'gene']
        results.append((gene, ';'.join(gene_ids) if gene_ids else 'NOT_FOUND'))
    except Exception as e:
        results.append((gene, f'ERROR:{e}'))
    time.sleep(0.5)

with open('data/reference/gene_symbol_to_id.tsv', 'w') as f:
    f.write("gene_symbol\tensembl_gene_id\n")
    for gene, gid in results:
        f.write(f"{gene}\t{gid}\n")

print(open('data/reference/gene_symbol_to_id.tsv').read())
