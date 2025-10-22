"""
GenomicsLab Pro - Automated Backend Test Script
------------------------------------------------
This script checks all 15 backend endpoints sequentially,
sends sample requests, validates the response, and prints
a color-coded summary report.

Run this while your backend (http://127.0.0.1:8000) is active.
"""

import requests
import time
import json
from colorama import Fore, Style, init
init(autoreset=True)

BASE_URL = "http://127.0.0.1:8000/sequence"

# --- Helper Function ---
def test_endpoint(name, path, payload):
    url = f"{BASE_URL}{path}"
    start = time.time()
    try:
        res = requests.post(url, json=payload, timeout=10)
        elapsed = round(time.time() - start, 2)
        if res.status_code == 200:
            print(Fore.GREEN + f"✅ {name:<30} → OK ({elapsed}s)")
            try:
                data = res.json()
                print(Fore.CYAN + f"   Sample Response: {json.dumps(data, indent=2)[:180]}...")
            except Exception:
                print(Fore.YELLOW + f"   ⚠️ Response not JSON-decodable")
        else:
            print(Fore.RED + f"❌ {name:<30} → HTTP {res.status_code}")
    except Exception as e:
        print(Fore.RED + f"❌ {name:<30} → FAILED ({e})")


if __name__ == "__main__":
    print(Style.BRIGHT + Fore.MAGENTA + "\n🚀 Running GenomicsLab Pro Backend Test Suite...\n")

    tests = [
        ("GC Content", "/gc_content", {"sequence": "ATGCGCGGATATAGCGC"}),
        ("Transcription", "/transcribe", {"sequence": "ATGCGTTA"}),
        ("Translation", "/translate", {"sequence": "AUGGCUUAA"}),
        ("Reverse Complement", "/reverse_complement", {"sequence": "ATGCGTTA"}),
        ("ORF Finder", "/orf_finder", {"sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"}),
        ("GC Sliding Window", "/gc_sliding_window", {"sequence": "ATGCGCGCGATATATATGCGCGCGTATATAT", "window_size": 10, "step_size": 5}),
        ("Codon Usage", "/codon_usage", {"sequence": "AUGGCUAUGGCUAUG"}),
        ("Restriction Sites", "/restriction_sites", {"sequence": "GAATTCGGATCCGCGGCCGCAAGCTT"}),
        ("Primer Design", "/primer_design", {"target_sequence": "ATGCGCGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT", "desired_tm": 60.0}),
        ("Motif Search", "/motif_search", {"sequence": "ATGCGCATGCATGCGC", "motif_pattern": "ATG"}),
        ("K-mer Analysis", "/kmer_analysis", {"sequence": "ATGCGATGCGATG", "k_length": 3}),
        ("Similarity Search", "/similarity_search", {"query_sequence": "ATGC", "target_sequence": "ATGCGCGCATGC", "alignment_type": "local"}),
        ("RNA Structure Prediction", "/rna_structure_prediction", {"sequence": "AUGGCUAUGCUAGCUAGCUGA"}),
        ("Signal Peptide Predictor", "/signal_peptide_predictor", {"sequence": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"}),
        ("Transmembrane Domain Finder", "/transmembrane_domain_finder", {"sequence": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEYGFQNALIVRYTRKPVSG"})
    ]

    for name, path, payload in tests:
        test_endpoint(name, path, payload)
        print()

    print(Style.BRIGHT + Fore.GREEN + "✅ All endpoints tested. Review output above for details.\n")
