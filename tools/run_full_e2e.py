#!/usr/bin/env python3
"""
Full end-to-end test driver for spades-rs v1.1.0-alpha.1.

Assembles every real-world dataset, runs QUAST where needed, then runs every
biological audit script against the resulting contigs/scaffolds and prints a
final scorecard and detailed report.

Usage:
    python3 tools/run_full_e2e.py [--skip-assemble]

    --skip-assemble  reuse previous assembly outputs (useful when re-running
                     audits without re-assembling)
"""

import sys
import os
import time
import subprocess
import shutil

ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN    = os.path.join(ROOT, "target", "release", "spades-rs")
EVAL   = os.path.join(ROOT, "target", "release", "spades-eval")
QUAST  = os.path.expanduser("~/.local/bin/quast.py")
DATA   = os.path.join(ROOT, "data")
OUT    = os.path.join(ROOT, "test_e2e_out")
TOOLS  = os.path.join(ROOT, "tools")

RESET  = "\033[0m"
BOLD   = "\033[1m"
GREEN  = "\033[32m"
RED    = "\033[31m"
YELLOW = "\033[33m"
CYAN   = "\033[36m"

# ---------------------------------------------------------------------------
# Assembly jobs
# ---------------------------------------------------------------------------
# Each entry: (name, cmd_args, output_subdir, ref_fasta)
ASSEMBLIES = [
    (
        "E. coli MG1655 (Illumina PE + Nanopore)",
        [BIN,
         "-1", f"{DATA}/ecoli_mg1655/SRR1770413_1.fastq.gz",
         "-2", f"{DATA}/ecoli_mg1655/SRR1770413_2.fastq.gz",
         "--nanopore", f"{DATA}/ecoli_mg1655/DRR242214_nanopore.fastq.gz",
         "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/ecoli"],
        f"{OUT}/ecoli",
        f"{DATA}/ecoli_mg1655/ref_NC_000913.3.fa",
    ),
    (
        "S. cerevisiae S288C (Illumina PE + ONT)",
        [BIN,
         "-1", f"{DATA}/yeast/yeast_R1.fq.gz",
         "-2", f"{DATA}/yeast/yeast_R2.fq.gz",
         "--nanopore", f"{DATA}/yeast/yeast_ont.fq.gz",
         "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/yeast"],
        f"{OUT}/yeast",
        f"{DATA}/yeast/s288c_ref.fa",
    ),
    (
        "S. pombe 972h- (Illumina PE)",
        [BIN,
         "-1", f"{DATA}/pombe/pombe_R1.fq.gz",
         "-2", f"{DATA}/pombe/pombe_R2.fq.gz",
         "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/pombe"],
        f"{OUT}/pombe",
        f"{DATA}/pombe/pombe_ref.fa",
    ),
    (
        "M. tuberculosis H37Rv (Illumina PE)",
        [BIN,
         "-1", f"{DATA}/tb/tb_R1.fq.gz",
         "-2", f"{DATA}/tb/tb_R2.fq.gz",
         "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/tb"],
        f"{OUT}/tb",
        f"{DATA}/tb/h37rv_ref.fa",
    ),
    (
        "P. aeruginosa PAO1 (Illumina PE)",
        [BIN,
         "-1", f"{DATA}/pa/pa_R1.fq.gz",
         "-2", f"{DATA}/pa/pa_R2.fq.gz",
         "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/pa"],
        f"{OUT}/pa",
        f"{DATA}/pa/pao1_ref.fa",
    ),
    (
        "P. falciparum 3D7 (Illumina PE)",
        [BIN,
         "-1", f"{DATA}/pf/pf_R1.fq.gz",
         "-2", f"{DATA}/pf/pf_R2.fq.gz",
         "-k", "33", "--multik", "21,33,55",
         "-o", f"{OUT}/pf"],
        f"{OUT}/pf",
        f"{DATA}/pf/pf3d7_ref.fa",
    ),
    (
        "ZymoBIOMICS mock community (meta, Illumina PE)",
        [BIN,
         "-1", f"{DATA}/zymo/zymo_R1.fq.gz",
         "-2", f"{DATA}/zymo/zymo_R2.fq.gz",
         "--meta", "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/zymo"],
        f"{OUT}/zymo",
        f"{DATA}/zymo/zymo_combined_ref.fa",
    ),
    (
        "K. pneumoniae BAA-2146 plasmids (--plasmid)",
        [BIN,
         "-1", f"{DATA}/plasmid_kpn/kpn_R1.fq.gz",
         "-2", f"{DATA}/plasmid_kpn/kpn_R2.fq.gz",
         "--plasmid", "-k", "55", "--multik", "33,55,77",
         "-o", f"{OUT}/plasmid_kpn"],
        f"{OUT}/plasmid_kpn",
        f"{DATA}/plasmid_kpn/kpn_baa2146_all_ref.fa",
    ),
    (
        "E. coli single-cell MDA (--sc)",
        [BIN,
         "-1", f"{DATA}/single_cell/sc_ecoli_R1.fq.gz",
         "-2", f"{DATA}/single_cell/sc_ecoli_R2.fq.gz",
         "--sc", "-k", "55", "--multik", "21,33,55",
         "-o", f"{OUT}/single_cell"],
        f"{OUT}/single_cell",
        f"{DATA}/single_cell/sc_ecoli_ref.fa",
    ),
    (
        "S. cerevisiae RNA-Seq transcriptome (--rna)",
        [BIN,
         "-1", f"{DATA}/rna/rna_R1.fq.gz",
         "-2", f"{DATA}/rna/rna_R2.fq.gz",
         "--rna", "-k", "33", "--multik", "21,33,55",
         "-o", f"{OUT}/rna"],
        f"{OUT}/rna",
        f"{DATA}/rna/scerevisiae_cdna.fa",
    ),
]

# ---------------------------------------------------------------------------
# QUAST alignment jobs needed by specific audit scripts
# ---------------------------------------------------------------------------
QUAST_JOBS = [
    (
        "QUAST: S. cerevisiae S288C",
        f"{OUT}/yeast/scaffolds.fasta",
        f"{DATA}/yeast/s288c_ref.fa",
        f"{OUT}/yeast/quast",
    ),
    (
        "QUAST: M. tuberculosis H37Rv",
        f"{OUT}/tb/scaffolds.fasta",
        f"{DATA}/tb/h37rv_ref.fa",
        f"{OUT}/tb/quast",
    ),
    (
        "QUAST: ZymoBIOMICS Contigs",
        f"{OUT}/zymo/contigs.fasta",
        f"{DATA}/zymo/zymo_combined_ref.fa",
        f"{OUT}/zymo/quast",
    ),
]

# ---------------------------------------------------------------------------
# Audit jobs  (name, script, extra_args_fn)
# ---------------------------------------------------------------------------
AUDITS = [
    (
        "E. coli rRNA operon synteny (7 operons)",
        "audit_rrna_synteny.py",
        lambda: [
            f"{OUT}/ecoli/scaffolds.fasta",
            f"{DATA}/ecoli_mg1655/ref_NC_000913.3.fa",
            f"{OUT}/ecoli/quast",
        ],
    ),
    (
        "S. cerevisiae chromosome synteny (16 chr)",
        "audit_yeast_chromosomes.py",
        lambda: [
            f"{OUT}/yeast/quast/contigs_reports/all_alignments_scaffolds.tsv",
        ],
    ),
    (
        "S. pombe loci (3 chromosomes)",
        "audit_pombe_loci.py",
        lambda: [
            f"{DATA}/pombe/pombe_ref.fa",
            f"{OUT}/pombe/scaffolds.fasta",
        ],
    ),
    (
        "M. tuberculosis AMR & PE/PPE genes",
        "audit_tb_genes.py",
        lambda: [
            f"{OUT}/tb/quast/contigs_reports/all_alignments_scaffolds.tsv",
        ],
    ),
    (
        "P. aeruginosa efflux & virulence genes",
        "audit_pa_genes.py",
        lambda: [
            f"{DATA}/pa/pao1_ref.fa",
            f"{DATA}/pa/pao1_annotations.gff",
            f"{OUT}/pa/scaffolds.fasta",
        ],
    ),
    (
        "P. falciparum 14-chromosome AT-rich genome",
        "audit_pf_genes.py",
        lambda: [
            f"{DATA}/pf/pf3d7_ref.fa",
            f"{DATA}/pf/pf3d7_annotations.gff",
            f"{OUT}/pf/scaffolds.fasta",
        ],
    ),
    (
        "ZymoBIOMICS per-species genome fraction & 16S/18S",
        "audit_zymo_community.py",
        lambda: [
            f"{DATA}/zymo/references/ZymoBIOMICS.STD.refseq.v2",
            f"{OUT}/zymo/scaffolds.fasta",
        ],
    ),
    (
        "ZymoBIOMICS chimera check",
        "check_zymo_chimeras.py",
        lambda: [
            f"{OUT}/zymo/quast/contigs_reports/all_alignments_contigs.tsv",
        ],
    ),
    (
        "K. pneumoniae plasmid circularity & AMR",
        "audit_plasmid_kpn.py",
        lambda: [
            f"{DATA}/plasmid_kpn/kpn_baa2146_all_ref.fa",
            f"{OUT}/plasmid_kpn/contigs.fasta",
            f"{OUT}/plasmid_kpn/plasmids.fasta",
        ],
    ),
    (
        "E. coli single-cell MDA coverage uniformity",
        "audit_single_cell_ecoli.py",
        lambda: [
            f"{DATA}/single_cell/sc_ecoli_ref.fa",
            f"{OUT}/single_cell/scaffolds.fasta",
        ],
    ),
    (
        "S. cerevisiae RNA-Seq transcript recovery",
        "audit_transcriptome_yeast.py",
        lambda: [
            f"{DATA}/rna/scerevisiae_cdna.fa",
            f"{OUT}/rna/scaffolds.fasta",
        ],
    ),
]


def run(label, cmd, log_path):
    """Run a subprocess, stream to log, return (returncode, elapsed)."""
    t0 = time.time()
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as lf:
        proc = subprocess.run(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=ROOT,
        )
    return proc.returncode, round(time.time() - t0, 1)


def print_header(title):
    print(f"\n{CYAN}{BOLD}{'=' * 70}{RESET}")
    print(f"{CYAN}{BOLD}  {title}{RESET}")
    print(f"{CYAN}{BOLD}{'=' * 70}{RESET}", flush=True)


def main():
    skip_assemble = "--skip-assemble" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    log_dir = os.path.join(OUT, "logs")
    os.makedirs(log_dir, exist_ok=True)

    results = []

    # -----------------------------------------------------------------------
    # Phase 1 — assemblies
    # -----------------------------------------------------------------------
    print_header("PHASE 1 — ASSEMBLIES")
    for name, cmd, out_dir, ref_fa in ASSEMBLIES:
        scaffolds = os.path.join(out_dir, "scaffolds.fasta")
        if skip_assemble and os.path.exists(scaffolds):
            print(f"  {YELLOW}[SKIP]{RESET} {name} (existing output reused)", flush=True)
            results.append({"phase": "assemble", "name": name, "passed": True,
                            "elapsed": 0, "skipped": True})
            continue

        # Check input files exist
        input_files = [a for a in cmd if a.endswith((".fq.gz", ".fastq.gz")) and not a.startswith("-")]
        missing = [f for f in input_files if not os.path.exists(f)]
        if missing:
            print(f"  {RED}[SKIP]{RESET} {name} — missing inputs: {', '.join(missing)}", flush=True)
            results.append({"phase": "assemble", "name": name, "passed": False,
                            "elapsed": 0, "skipped": True,
                            "reason": f"missing: {', '.join(os.path.basename(m) for m in missing)}"})
            continue

        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        safe = name.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")
        log = os.path.join(log_dir, f"assemble_{safe}.log")
        print(f"  {BOLD}▶{RESET} {name} ...", flush=True)
        rc, elapsed = run(name, cmd, log)
        passed = rc == 0 and os.path.exists(scaffolds)
        
        # Run spades-eval to append assembly metrics to log
        if passed and os.path.exists(EVAL) and ref_fa and os.path.exists(ref_fa):
            with open(log, "a") as lf:
                lf.write("\n--- SPADES-EVAL SUMMARY ---\n")
                subprocess.run([EVAL, scaffolds, "-r", ref_fa], stdout=lf, stderr=subprocess.STDOUT)

        tag = f"{GREEN}PASS{RESET}" if passed else f"{RED}FAIL{RESET}"
        print(f"  [{tag}] {name:<55} {elapsed}s", flush=True)
        if not passed:
            print(f"         ↳ log: {log}", flush=True)
        results.append({"phase": "assemble", "name": name, "passed": passed,
                        "elapsed": elapsed, "log": log})

    # -----------------------------------------------------------------------
    # Phase 1.5 — QUAST alignments for audit scripts
    # -----------------------------------------------------------------------
    if os.path.exists(QUAST):
        print_header("PHASE 1.5 — QUAST REFERENCE ALIGNMENTS")
        for qname, c_fa, r_fa, q_out in QUAST_JOBS:
            if not os.path.exists(c_fa) or not os.path.exists(r_fa):
                continue
            safe = qname.replace(" ", "_").replace(":", "")
            log = os.path.join(log_dir, f"quast_{safe}.log")
            print(f"  {BOLD}▶{RESET} {qname} ...", flush=True)
            cmd = ["python3", QUAST, c_fa, "-r", r_fa, "-o", q_out, "--threads", "22", "--min-contig", "200"]
            rc, elapsed = run(qname, cmd, log)
            tag = f"{GREEN}DONE{RESET}" if rc == 0 else f"{YELLOW}WARN{RESET}"
            print(f"  [{tag}] {qname:<55} {elapsed}s", flush=True)

    # -----------------------------------------------------------------------
    # Phase 2 — biological audits
    # -----------------------------------------------------------------------
    print_header("PHASE 2 — BIOLOGICAL AUDITS")
    for name, script, args_fn in AUDITS:
        args = args_fn()
        # Check input paths exist (skip output dirs like quast if first arg is input)
        check_paths = args[:2] if script == "audit_rrna_synteny.py" else args
        missing = [a for a in check_paths if not os.path.exists(a)]
        if missing:
            print(f"  {YELLOW}[SKIP]{RESET} {name} — missing: {', '.join(os.path.basename(m) for m in missing)}", flush=True)
            results.append({"phase": "audit", "name": name, "passed": False,
                            "elapsed": 0, "skipped": True,
                            "reason": f"missing: {', '.join(os.path.basename(m) for m in missing)}"})
            continue

        cmd = ["python3", os.path.join(TOOLS, script)] + args
        safe = name.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "").replace("/", "_")
        log = os.path.join(log_dir, f"audit_{safe}.log")
        print(f"  {BOLD}▶{RESET} {name} ...", flush=True)
        rc, elapsed = run(name, cmd, log)
        passed = rc == 0
        tag = f"{GREEN}PASS{RESET}" if passed else f"{RED}FAIL{RESET}"
        print(f"  [{tag}] {name:<55} {elapsed}s", flush=True)
        if not passed:
            print(f"         ↳ log: {log}", flush=True)
        results.append({"phase": "audit", "name": name, "passed": passed,
                        "elapsed": elapsed, "log": log})

    # -----------------------------------------------------------------------
    # Scorecard
    # -----------------------------------------------------------------------
    print(f"\n{BOLD}{'=' * 70}{RESET}")
    print(f"{BOLD}  END-TO-END SCORECARD — spades-rs v1.1.0-alpha.1{RESET}")
    print(f"{BOLD}{'=' * 70}{RESET}")

    sections = [
        ("Assemblies",  [r for r in results if r["phase"] == "assemble"]),
        ("Bio audits",  [r for r in results if r["phase"] == "audit"]),
    ]

    grand_pass = grand_fail = 0
    for section, items in sections:
        if not items:
            continue
        print(f"\n  {BOLD}{section}{RESET}")
        for r in items:
            if r.get("skipped") and r["passed"]:
                tag = f"{YELLOW}SKIP{RESET}"
            elif r["passed"]:
                tag = f"{GREEN}PASS{RESET}"
                grand_pass += 1
            else:
                tag = f"{RED}FAIL{RESET}"
                grand_fail += 1
            label = r["name"][:60]
            elapsed = r.get("elapsed", 0)
            print(f"    [{tag}] {label:<62} {elapsed}s")
            if not r["passed"] and "reason" in r:
                print(f"           ↳ {r['reason']}")

    total = grand_pass + grand_fail
    print(f"\n{BOLD}{'─' * 70}{RESET}")
    print(f"  {BOLD}Passed:{RESET} {GREEN}{grand_pass}{RESET} / {total}   "
          f"{BOLD}Failed:{RESET} {RED}{grand_fail}{RESET}")
    print(f"{BOLD}{'=' * 70}{RESET}\n", flush=True)

    sys.exit(0 if grand_fail == 0 else 1)


if __name__ == "__main__":
    main()
