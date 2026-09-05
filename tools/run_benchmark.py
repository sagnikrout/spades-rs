#!/usr/bin/env python3
"""
Automated Benchmarking & Profiling Harness for Assemblers.
Executes an assembler pipeline under GNU /usr/bin/time -v profiler,
extracts Peak RSS and Timing metrics, runs evaluation, and logs results.
"""

import sys
import os
import subprocess
import re
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def parse_time_output(stderr_text):
    metrics = {}
    for line in stderr_text.splitlines():
        line = line.strip()
        if "User time (seconds):" in line:
            metrics["user_time_s"] = float(line.split(":")[-1].strip())
        elif "System time (seconds):" in line:
            metrics["sys_time_s"] = float(line.split(":")[-1].strip())
        elif "Percent of CPU this job got:" in line:
            metrics["cpu_percent"] = line.split(":")[-1].strip()
        elif "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
            raw_t = line.split("):")[-1].strip()
            parts = raw_t.split(":")
            if len(parts) == 3:
                secs = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                secs = float(parts[0]) * 60 + float(parts[1])
            else:
                secs = float(parts[0])
            metrics["wall_time_s"] = round(secs, 3)
            metrics["wall_time_str"] = raw_t
        elif "Maximum resident set size (kbytes):" in line:
            rss_kb = int(line.split(":")[-1].strip())
            metrics["peak_rss_kb"] = rss_kb
            metrics["peak_rss_mb"] = round(rss_kb / 1024.0, 2)
    return metrics

def run_benchmarked_command(name, command_str, output_fasta, ref_fasta=None):
    print("=" * 60)
    print(f"BENCHMARK RUN: {name}")
    print(f"Command: {command_str}")
    print("=" * 60)
    
    time_log = f"/tmp/bench_time_{os.getpid()}.log"
    wrapped_cmd = f"/usr/bin/time -v -o {time_log} {command_str}"
    
    t0 = time.time()
    res = subprocess.run(wrapped_cmd, shell=True)
    t1 = time.time()
    
    time_metrics = {}
    if os.path.exists(time_log):
        with open(time_log) as f:
            time_metrics = parse_time_output(f.read())
    else:
        time_metrics["wall_time_s"] = round(t1 - t0, 3)
        
    print("\n--- PROFILER RESOURCE USAGE ---")
    print(f"  Wall Time:       {time_metrics.get('wall_time_s', 'N/A')} s")
    print(f"  User CPU Time:   {time_metrics.get('user_time_s', 'N/A')} s")
    print(f"  Peak RSS:        {time_metrics.get('peak_rss_mb', 'N/A')} MB")
    if "cpu_percent" in time_metrics:
        print(f"  CPU Allocation:  {time_metrics['cpu_percent']}")
        
    # Evaluate assembly quality
    assembly_metrics = {}
    if os.path.exists(output_fasta):
        from eval_assembly import evaluate_assembly
        assembly_metrics = evaluate_assembly(output_fasta, ref_fasta)
        print("\n--- ASSEMBLY QUALITY ---")
        for k, v in assembly_metrics.items():
            print(f"  {k:<35}: {v}")
    else:
        print(f"WARNING: Output file {output_fasta} not found!")
        
    result = {
        "name": name,
        "exit_code": res.returncode,
        "profiler": time_metrics,
        "assembly": assembly_metrics,
    }
    return result

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run benchmark and profile assembly")
    parser.add_argument("--name", required=True, help="Benchmark run name")
    parser.add_argument("--cmd", required=True, help="Command to run")
    parser.add_argument("--output", required=True, help="Output contigs.fasta path")
    parser.add_argument("--ref", help="Optional reference FASTA")
    parser.add_argument("--save_json", help="Save metrics as JSON")
    
    args = parser.parse_args()
    res = run_benchmarked_command(args.name, args.cmd, args.output, args.ref)
    
    if args.save_json:
        if os.path.dirname(args.save_json):
            os.makedirs(os.path.dirname(args.save_json), exist_ok=True)
        with open(args.save_json, "w") as f:
            json.dump(res, f, indent=2)
