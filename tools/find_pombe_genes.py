import os
import sys

gff_path = sys.argv[1] if len(sys.argv) > 1 else 'data/pombe/pombe_annotations.gff'
if not os.path.exists(gff_path):
    print(f"Error: Annotations file '{gff_path}' not found.")
    sys.exit(1)

with open(gff_path) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) >= 9 and parts[0] == 'NC_003423.3':
            s = int(parts[3])
            e = int(parts[4])
            if 2115000 <= s <= 2126500:
                print(f"{parts[2]:12s} | {s:8d} - {e:8d} | {parts[8][:100]}")
