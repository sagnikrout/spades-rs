from audit_pombe_loci import load_fasta

ref = load_fasta('data/pombe/pombe_ref.fa')
scaffolds = load_fasta('output/bench_pombe/scaffolds.fasta')

cdc2 = ref['NC_003423.3'][1500208-1:1502095].upper()
s3 = None
for k, v in scaffolds.items():
    if 'scaffold_3_' in k:
        s3 = v.upper()
        break

start = 93574
s3_sub = s3[start : start + len(cdc2)]

diffs = sum(1 for a, b in zip(cdc2, s3_sub) if a != b)
print(f"CDC2 uppercase comparison over {len(cdc2)} bp:")
print(f"Mismatches: {diffs} / {len(cdc2)} (Identity: {(len(cdc2)-diffs)/len(cdc2)*100:.2f}%)")
