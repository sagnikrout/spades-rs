with open('data/pombe/pombe_annotations.gff') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) >= 9 and parts[0] == 'NC_003423.3':
            s = int(parts[3])
            e = int(parts[4])
            if 2115000 <= s <= 2126500:
                print(f"{parts[2]:12s} | {s:8d} - {e:8d} | {parts[8][:100]}")
