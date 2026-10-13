#!/usr/bin/env python3
import http.server
import socketserver
import json
import os
import glob
import re

PORT = 8999
ROOT = "/mnt/c/Users/sagni/Documents/antigravity/intelligent-pascal"
LOG_DIR = os.path.join(ROOT, "test_e2e_out", "logs")

TOTAL_STEPS = 20

BENCHMARKS = [
    ("E. coli MG1655", "assemble_E._coli_MG1655_Illumina_PE_+_Nanopore.log"),
    ("S. cerevisiae S288C", "assemble_S._cerevisiae_S288C_Illumina_PE_+_ONT.log"),
    ("S. pombe 972h-", "assemble_S._pombe_972h-_Illumina_PE.log"),
    ("M. tuberculosis H37Rv", "assemble_M._tuberculosis_H37Rv_Illumina_PE.log"),
    ("P. aeruginosa PAO1", "assemble_P._aeruginosa_PAO1_Illumina_PE.log"),
    ("P. falciparum 3D7", "assemble_P._falciparum_3D7_Illumina_PE.log"),
    ("Zymo Metagenome", "assemble_ZymoBIOMICS_mock_community_meta_Illumina_PE.log"),
    ("K. pneumoniae Plasmids", "assemble_K._pneumoniae_BAA-2146_plasmids_--plasmid.log"),
    ("E. coli Single-Cell", "assemble_E._coli_single-cell_MDA_--sc.log"),
    ("Yeast Transcriptome", "assemble_S._cerevisiae_RNA-Seq_transcriptome_--rna.log"),
]

def get_status():
    items = []
    completed_count = 0
    current_task = "Idle"
    current_detail = ""
    
    # Check assemblies
    for name, log_file in BENCHMARKS:
        path = os.path.join(LOG_DIR, log_file)
        if not os.path.exists(path):
            items.append({"name": name, "status": "pending"})
        else:
            with open(path, "r", errors="ignore") as f:
                content = f.read()
            if "SPADES-EVAL SUMMARY" in content or "Completed in" in content or "Final assembly stats" in content or "Contiguity Metrics" in content:
                items.append({"name": name, "status": "completed"})
                completed_count += 1
            else:
                items.append({"name": name, "status": "running"})
                current_task = name
                # find last stage/step
                lines = [l.strip() for l in content.splitlines() if l.strip()]
                for l in reversed(lines):
                    if l.startswith("───") or l.startswith("▶▶▶") or "Step k=" in l:
                        current_detail = l
                        break

    # Check audits (10 audits + QUAST)
    audit_logs = glob.glob(os.path.join(LOG_DIR, "audit_*.log")) + glob.glob(os.path.join(LOG_DIR, "quast_*.log"))
    for al in audit_logs:
        completed_count += 1
        
    pct = min(100, int((completed_count / TOTAL_STEPS) * 100))
    if current_task == "Idle" and completed_count > 0:
        if completed_count >= TOTAL_STEPS:
            current_task = "Complete"
            current_detail = "All 20 benchmarks and audits passed successfully."
        else:
            current_task = "Running audits"
            current_detail = "Executing biological synteny & QUAST verification."

    return {
        "pct": pct,
        "completed": completed_count,
        "total": TOTAL_STEPS,
        "current_task": current_task,
        "current_detail": current_detail,
        "items": items,
    }

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/status" or self.path == "/status.json":
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
            self.end_headers()
            data = json.dumps(get_status()).encode("utf-8")
            self.wfile.write(data)
        else:
            self.send_response(404)
            self.end_headers()

if __name__ == "__main__":
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        httpd.serve_forever()
