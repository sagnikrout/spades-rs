import glob
import os

ref_dir = 'data/zymo/references/ZymoBIOMICS.STD.refseq.v2/Genomes/'
files = sorted(glob.glob(os.path.join(ref_dir, '*.fasta')))

out_path = 'data/zymo/zymo_combined_ref.fa'
total_bp = 0
species_counts = {}

with open(out_path, 'w') as out_f:
    for f in files:
        bname = os.path.basename(f)
        sp_name = bname.replace('_complete_genome.fasta', '').replace('_draft_genome.fasta', '')
        sp_bp = 0
        chr_idx = 1
        with open(f, 'r') as in_f:
            for line in in_f:
                if line.startswith('>'):
                    orig_h = line[1:].strip()
                    out_f.write(f">{sp_name}_contig_{chr_idx} {orig_h}\n")
                    chr_idx += 1
                else:
                    seq = line.strip()
                    sp_bp += len(seq)
                    out_f.write(f"{seq}\n")
        species_counts[sp_name] = sp_bp
        total_bp += sp_bp
        print(f"  {sp_name:25s}: {sp_bp:10,d} bp")

print("=" * 60)
print(f"Combined reference saved to: {out_path}")
print(f"Total reference genome size: {total_bp:,} bp across {len(species_counts)} species.")
