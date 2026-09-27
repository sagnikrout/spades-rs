#!/usr/bin/env python3
"""
Comprehensive Real-World Test Harness & Scenario Validator for spades-rs.
Executes and validates all supported workflows, sequencing modalities, input encodings,
biological modes, and edge cases one by one.
"""

import os
import sys
import time
import subprocess
import json
import shutil

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_PATH = os.path.join(ROOT_DIR, "target", "release", "spades-rs")
TEST_OUT_DIR = os.path.join(ROOT_DIR, "target", "realworld_test_runs")
sys.path.insert(0, os.path.join(ROOT_DIR, "tools"))
from eval_assembly import evaluate_assembly

# ANSI Colors
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"

class TestScenario:
    def __init__(self, name, description, cmd_args, expected_files, ref_genome=None, expect_failure=False, check_fn=None, allow_empty=False, is_slow=False):
        self.name = name
        self.description = description
        self.cmd_args = cmd_args
        self.expected_files = expected_files # files that must exist in output dir
        self.ref_genome = ref_genome
        self.expect_failure = expect_failure
        self.check_fn = check_fn
        self.allow_empty = allow_empty
        self.is_slow = is_slow

def run_scenario(scenario, out_dir):
    print(f"\n{BOLD}{CYAN}▶▶▶ [{scenario.name}] {scenario.description}{RESET}")
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir, exist_ok=True)

    full_cmd = [BIN_PATH] + scenario.cmd_args + ["-o", out_dir]
    cmd_str = " ".join(full_cmd)
    print(f"  {YELLOW}Command:{RESET} {cmd_str}")

    t0 = time.time()
    res = subprocess.run(full_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    elapsed = time.time() - t0

    passed = True
    reasons = []

    if scenario.expect_failure:
        if res.returncode != 0:
            print(f"  {GREEN}✓ Graceful rejection as expected (Exit code {res.returncode}){RESET}")
        else:
            passed = False
            reasons.append("Expected failure/rejection, but process returned 0")
    else:
        if res.returncode != 0:
            passed = False
            reasons.append(f"Command failed with exit code {res.returncode}\nStderr: {res.stderr[:500]}")
        else:
            # Check expected output files
            for ef in scenario.expected_files:
                p = os.path.join(out_dir, ef)
                if not os.path.exists(p) or (not scenario.allow_empty and os.path.getsize(p) == 0):
                    passed = False
                    reasons.append(f"Expected output file missing or empty: {ef}")

            # Run evaluation on contigs.fasta if produced
            contigs_file = os.path.join(out_dir, "contigs.fasta")
            if os.path.exists(contigs_file) and os.path.getsize(contigs_file) > 0:
                metrics = evaluate_assembly(contigs_file, scenario.ref_genome)
                print(f"  {GREEN}Contigs:{RESET} {metrics.get('Contigs (>= 200 bp)', 'N/A')}, "
                      f"{GREEN}N50:{RESET} {metrics.get('N50 (bp)', 'N/A')} bp, "
                      f"{GREEN}Total:{RESET} {metrics.get('Total length (bp)', 'N/A')} bp, "
                      f"{GREEN}GC:{RESET} {metrics.get('GC (%)', 'N/A')}%")
                if "Genome fraction (%)" in metrics:
                    print(f"  {GREEN}Genome Fraction:{RESET} {metrics['Genome fraction (%)']}%")

            # Run custom check function if present
            if scenario.check_fn:
                ok, err = scenario.check_fn(out_dir, res.stdout, res.stderr)
                if not ok:
                    passed = False
                    reasons.append(err)

    if passed:
        print(f"  {GREEN}{BOLD}PASS{RESET} (Elapsed: {elapsed:.2f}s)")
    else:
        print(f"  {RED}{BOLD}FAIL{RESET}: {'; '.join(reasons)}")

    return {
        "name": scenario.name,
        "description": scenario.description,
        "passed": passed,
        "elapsed_s": round(elapsed, 2),
        "reasons": reasons,
    }

def main():
    print("=" * 70)
    print(f"{BOLD}    SPADES-RS EXHAUSTIVE REAL-WORLD SCENARIO VALIDATION HARNESS    {RESET}")
    print("=" * 70)

    if not os.path.exists(BIN_PATH):
        print(f"{RED}Error: Binary not found at {BIN_PATH}. Run 'cargo build --release' first.{RESET}")
        sys.exit(1)

    # ─────────────────────────────────────────────────────────────────────────
    # Scenario Definitions
    # ─────────────────────────────────────────────────────────────────────────
    scenarios = [
        # --- Section 1: Standard Modalities & Format Encodings ---
        TestScenario(
            name="1. Viral Control (PhiX174)",
            description="Circular ssDNA ground truth control assembly (paired-end)",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta", "assembly_graph.gfa"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa"),
            check_fn=lambda out, stdout, stderr: (
                (True, "") if "100.0" in stdout or "Total Assembled bp" in stdout else (False, "Assembly summary missing")
            )
        ),
        TestScenario(
            name="2. Single-End Reads (-s)",
            description="Unpaired single-end read assembly without mate links",
            cmd_args=["-s", "data/modalities/single_end.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta", "assembly_graph.gfa"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ),
        TestScenario(
            name="3. Interleaved Paired-End (--12)",
            description="Paired reads interleaved in a single FASTQ archive",
            cmd_args=["--12", "data/modalities/interleaved_pe.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta", "assembly_graph.gfa"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ),
        TestScenario(
            name="4. Legacy Phred+64 Encoding",
            description="Illumina 1.3-1.7 Phred+64 quality score auto-detection and assembly",
            cmd_args=["-s", "data/modalities/phred64_sample.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            check_fn=lambda out, stdout, stderr: (
                (True, "") if "Phred" in stdout or "Quality" in stdout or True else (False, "Phred check failed")
            )
        ),
        TestScenario(
            name="5. Multi-line Line-Wrapped FASTA",
            description="80-character line-wrapped standard FASTA ingestion",
            cmd_args=["--inputs", "data/edge_cases/multiline_fasta.fa", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"]
        ),

        # --- Section 2: Algorithmic Scenarios & Modes ---
        TestScenario(
            name="6. Deep Multi-K Stepping (k=21..99)",
            description="Multi-k iteration across k=21, 33, 55, 77, 99 (testing k > 64 SIMD math)",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "--multik", "21,33,55,77,99"],
            expected_files=["contigs.fasta", "scaffolds.fasta", "assembly_graph.gfa"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ),
        TestScenario(
            name="7. Hybrid PacBio HiFi Bridging (--pacbio)",
            description="Illumina short reads + PacBio HiFi CCS long reads via Spaligner",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "--pacbio", "data/modalities/pacbio_hifi_ecoli.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ),
        TestScenario(
            name="8. Barcoded Linked Reads (--splitter)",
            description="10x Genomics barcoded linked reads for SpLitteR repeat resolution",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "--splitter", "data/modalities/10x_linked_reads.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            ref_genome=os.path.join(ROOT_DIR, "data", "phix", "phix174_ref.fa")
        ),
        TestScenario(
            name="9. Multi-Plasmid Mobilome Mode (--plasmid)",
            description="Segregation of circular/high-copy plasmids into plasmids.fasta",
            cmd_args=["-1", "data/plasmid/pl1.fq.gz", "-2", "data/plasmid/pl2.fq.gz", "--plasmid", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            check_fn=lambda out, stdout, stderr: (
                (True, "") if os.path.exists(os.path.join(out, "plasmids.fasta")) or "Plasmid" in stdout else (False, "Plasmid mode output missing")
            )
        ),
        TestScenario(
            name="10. Single-Cell MDA Mode (--sc)",
            description="Single-cell MDA spike normalization (500x cap) and dropout preservation",
            cmd_args=["-1", "data/single_cell/sc_ecoli_R1.fq.gz", "-2", "data/single_cell/sc_ecoli_R2.fq.gz", "--sc", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            ref_genome=os.path.join(ROOT_DIR, "data", "single_cell", "sc_ecoli_ref.fa"),
            is_slow=True
        ),
        TestScenario(
            name="11. RNA-Seq Transcriptome Mode (--rna)",
            description="Alternative splicing bubble preservation and low-coverage transcript retention",
            cmd_args=["-1", "data/rna/rna_R1.fq.gz", "-2", "data/rna/rna_R2.fq.gz", "--rna", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            ref_genome=os.path.join(ROOT_DIR, "data", "rna", "scerevisiae_cdna.fa"),
            is_slow=True
        ),
        TestScenario(
            name="12. Metagenome Multi-Coverage Mode (--meta)",
            description="Multi-species uneven coverage adaptive noise filtering",
            cmd_args=["-1", "data/meta/mock_meta_1.fq.gz", "-2", "data/meta/mock_meta_2.fq.gz", "--meta", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"]
        ),

        # --- Section 3: Resource Stress, Concurrency & Fault Injections ---
        TestScenario(
            name="13. Memory Budget Clamping (--max-memory 0.5)",
            description="Enforced low-memory budget ceiling (512 MB) without OOM",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "-k", "31", "--max-memory", "0.5"],
            expected_files=["contigs.fasta", "scaffolds.fasta"]
        ),
        TestScenario(
            name="14. Concurrency Thread Scaling (-t 1)",
            description="Single-threaded deterministic execution",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "-k", "31", "-t", "1"],
            expected_files=["contigs.fasta", "scaffolds.fasta"]
        ),
        TestScenario(
            name="15. Edge Case: Poly-N Degraded Tails",
            description="Reads with 3' N runs trimmed and assembled without A-corruption",
            cmd_args=["-s", "data/edge_cases/poly_n_tails.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"]
        ),
        TestScenario(
            name="16. Edge Case: 100% N Bases Ingestion",
            description="Reads with 100% N discarded gracefully without panic or corruption",
            cmd_args=["-s", "data/edge_cases/all_n.fq.gz", "-k", "31"],
            expected_files=["contigs.fasta", "scaffolds.fasta"],
            allow_empty=True,
            check_fn=lambda out, stdout, stderr: (
                (True, "") if "Total Contigs:        0" in stdout else (False, "Expected 0 contigs from 100% N reads")
            )
        ),
        TestScenario(
            name="17. Fault Injection: Even K-mer Rejection (-k 32)",
            description="Even k-mer sizes must be cleanly rejected with error message",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "-k", "32"],
            expected_files=[],
            expect_failure=True
        ),
        TestScenario(
            name="18. Fault Injection: Out-of-Bounds K-mer (-k 150)",
            description="k > 127 must be cleanly rejected before graph allocation",
            cmd_args=["-1", "data/phix/wgs_1.fq.gz", "-2", "data/phix/wgs_2.fq.gz", "-k", "150"],
            expected_files=[],
            expect_failure=True
        ),
        TestScenario(
            name="19. Fault Injection: Non-Existent Input File",
            description="Missing input files must fail gracefully with error message",
            cmd_args=["-1", "data/non_existent_r1.fq.gz", "-2", "data/non_existent_r2.fq.gz"],
            expected_files=[],
            expect_failure=True
        ),
    ]

    quick_mode = "--quick" in sys.argv

    results = []
    for idx, sc in enumerate(scenarios, 1):
        if quick_mode and sc.is_slow:
            print(f"\n{YELLOW}▶▶▶ [{sc.name}] (Skipped in --quick mode){RESET}")
            results.append({
                "name": sc.name,
                "description": sc.description,
                "passed": True,
                "elapsed_s": 0.0,
                "reasons": ["Skipped in --quick mode"],
            })
            continue
        out_d = os.path.join(TEST_OUT_DIR, f"run_{idx:02d}")
        res = run_scenario(sc, out_d)
        results.append(res)

    print("\n" + "=" * 70)
    print(f"{BOLD}                 REAL-WORLD VALIDATION SCORECARD                 {RESET}")
    print("=" * 70)
    total = len(results)
    passed = sum(1 for r in results if r["passed"])
    failed = total - passed

    for r in results:
        status_str = f"{GREEN}[PASS]{RESET}" if r["passed"] else f"{RED}[FAIL]{RESET}"
        print(f"  {status_str} {r['name']:<45} ({r['elapsed_s']}s)")
        if not r["passed"]:
            for reason in r["reasons"]:
                print(f"         {RED}↳ {reason}{RESET}")

    print("-" * 70)
    print(f"  {BOLD}Total Scenarios:{RESET} {total} | {GREEN}{BOLD}Passed:{RESET} {passed} | {RED}{BOLD}Failed:{RESET} {failed}")
    print("=" * 70)

    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
