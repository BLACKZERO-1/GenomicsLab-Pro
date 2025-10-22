# backend/sequence_processing/core.py
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, EcoRI, BamHI, HindIII, NotI
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np
import re
from Bio import pairwise2
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Helper function (single definition)
def get_codon_usage_data(seq: str) -> dict:
    codon_counts = {}
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
    return codon_counts

# ----------------- 1. GC Content Calculator -----------------
def calculate_gc_content(dna_sequence: str) -> float:
    """Calculates the GC content percentage of a DNA/RNA sequence."""
    try:
        sequence = Seq(dna_sequence.upper())
        gc_frac = gc_fraction(sequence, ambiguous='remove')
        return round(gc_frac * 100.0, 2)
    except Exception:
        return 0.0

# ----------------- 2. DNA RNA Transcription -----------------
def transcribe_dna_to_rna(dna_sequence: str) -> str:
    """Transcribes a DNA sequence to its corresponding RNA sequence (T -> U)."""
    try:
        dna_seq_obj = Seq(dna_sequence.upper())
        rna_seq_obj = dna_seq_obj.transcribe()
        return str(rna_seq_obj)
    except Exception:
        return "ERROR: Invalid DNA sequence for transcription."

# ----------------- 3. RNA Protein Translation -----------------
def translate_rna_to_protein(rna_sequence: str) -> str:
    """
    Translates an RNA sequence to its corresponding amino acid protein sequence.
    """
    try:
        rna_seq_obj = Seq(rna_sequence.upper())
        protein_seq_obj = rna_seq_obj.translate(table=1, stop_symbol="*", to_stop=False)
        return str(protein_seq_obj)
    except Exception:
        return "ERROR: Invalid RNA sequence or frame for translation."

# ----------------- 4. Reverse Complement -----------------
def get_reverse_complement(dna_sequence: str) -> str:
    """Calculates the reverse complement of a DNA sequence."""
    try:
        seq_obj = Seq(dna_sequence.upper())
        rev_comp_obj = seq_obj.reverse_complement()
        return str(rev_comp_obj)
    except Exception:
        return "ERROR: Invalid DNA sequence for reverse complement."

# ----------------- 5. ORF Finder -----------------
def find_orfs(dna_sequence: str, min_peptide_length: int = 20) -> dict:
    """Finds potential Open Reading Frames (ORFs) in all six frames."""
    try:
        seq_obj = Seq(dna_sequence.upper())
        orfs = {}

        for frame in range(3):
            # Forward Frames
            peptide = str(seq_obj[frame:].translate(to_stop=False))
            potential_peptides = [p for p in peptide.split('*') if len(p) >= min_peptide_length]
            if potential_peptides:
                orfs[f"Frame +{frame + 1}"] = potential_peptides

            # Reverse Frames
            rev_comp_obj = seq_obj.reverse_complement()
            peptide_rev = str(rev_comp_obj[frame:].translate(to_stop=False))
            potential_peptides_rev = [p for p in peptide_rev.split('*') if len(p) >= min_peptide_length]
            if potential_peptides_rev:
                orfs[f"Frame -{frame + 1}"] = potential_peptides_rev

        return orfs

    except Exception:
        return {"ERROR": "Invalid DNA sequence or frame for ORF analysis."}

# ----------------- 6. GC Sliding Window -----------------
def gc_sliding_window(dna_sequence: str, window: int, step: int) -> list[float]:
    """Calculates the GC content over a sequence using a sliding window."""
    if window <= 0 or step <= 0 or window > len(dna_sequence):
        return []

    gc_list = []
    seq_upper = dna_sequence.upper()
    n = len(seq_upper)

    for i in range(0, n - window + 1, step):
        window_seq = seq_upper[i:i + window]
        gc_perc = calculate_gc_content(str(window_seq))
        gc_list.append(gc_perc)

    return gc_list

# ----------------- 7. Codon Usage Analysis -----------------
def analyze_codon_usage(rna_sequence: str) -> dict:
    """
    Calculates the usage frequency (per 1000 codons) and total count
    for each codon in a coding sequence.
    """
    try:
        rna_seq_upper = rna_sequence.upper()
        if len(rna_seq_upper) % 3 != 0:
            return {"ERROR": "Sequence length must be a multiple of 3 for codon analysis."}

        total_codons = len(rna_seq_upper) // 3
        if total_codons == 0:
            return {"ERROR": "Sequence is too short or invalid."}

        codon_counts = get_codon_usage_data(rna_seq_upper)
        usage_data = {}

        for codon, count in codon_counts.items():
            frequency = round((count / total_codons) * 1000, 2)
            usage_data[codon] = {
                "count": count,
                "frequency_per_1000": frequency
            }

        return {
            "total_codons": total_codons,
            "codon_usage": usage_data
        }

    except Exception:
        return {"ERROR": "Invalid RNA sequence for codon usage analysis."}

# ----------------- 8. Restriction Site Finder -----------------
def find_restriction_sites(dna_sequence: str) -> dict:
    """Finds cutting sites for a standard batch of restriction enzymes."""
    try:
        cleaned_sequence = "".join(c for c in dna_sequence.upper() if c in 'ATGC')
        if not cleaned_sequence:
            raise ValueError("Input sequence contained no valid DNA bases (A, T, G, C).")

        seq_obj = Seq(cleaned_sequence)

        batch = RestrictionBatch([EcoRI, BamHI, HindIII, NotI])

        analysis = batch.search(seq_obj)

        result = {}
        for enzyme, sites in analysis.items():
            if sites:
                result[enzyme.name] = [s + 1 for s in sites]

        if not result:
            return {"ERROR": "No cut sites found for tested enzymes."}

        return result

    except ValueError as e:
        return {"ERROR": f"Validation Error: {str(e)}"}
    except Exception as e:
        # keep error message for debugging but return a safe message to caller
        print(f"Restriction Site Analysis Failed: {e}")
        return {"ERROR": "Unexpected error during restriction analysis. Check sequence format."}

# ----------------- 9. Primer Design Tool -----------------
def design_pcr_primers(target_sequence: str, desired_tm: float = 60.0) -> dict:
    """Performs basic validation and estimates properties for a theoretical primer pair."""
    if len(target_sequence) < 50:
        return {"ERROR": "Target sequence must be at least 50 bases for meaningful primer design."}

    forward_seq = target_sequence[:18].upper()
    reverse_target = target_sequence[-18:].upper()
    reverse_seq = str(Seq(reverse_target).reverse_complement())

    try:
        tm_f = round(mt.Tm_NN(forward_seq), 2)
        tm_r = round(mt.Tm_NN(reverse_seq), 2)
        gc_f = calculate_gc_content(forward_seq)
        gc_r = calculate_gc_content(reverse_seq)

        tm_diff = round(abs(tm_f - tm_r), 2)

        return {
            "Forward_Primer": {
                "sequence": forward_seq,
                "Tm_NN": tm_f,
                "GC_Content": gc_f
            },
            "Reverse_Primer": {
                "sequence": reverse_seq,
                "Tm_NN": tm_r,
                "GC_Content": gc_r
            },
            "Summary": {
                "Length": 18,
                "Target_Tm": desired_tm,
                "Tm_Difference": tm_diff
            }
        }
    except Exception:
        return {"ERROR": "Calculation failed. Check sequence base integrity."}

# ----------------- 10. Motif Search -----------------
def find_sequence_motifs(sequence: str, motif_pattern: str, is_dna: bool = True) -> dict:
    """
    Searches for a user-defined pattern (motif) in a sequence.
    Motif pattern should be a standard IUPAC code or regular expression.
    """
    try:
        seq_upper = sequence.upper()
        if not seq_upper:
            return {"ERROR": "Input sequence is empty."}

        # 1. Handle common ambiguous IUPAC codes for DNA/RNA
        pattern = motif_pattern.upper()
        pattern = pattern.replace('W', '[AT]').replace('S', '[GC]').replace('R', '[AG]').replace('Y', '[CT]').replace('K', '[GT]').replace('M', '[AC]')

        # 2. Prepare search sequence according to DNA/RNA
        if is_dna:
            search_seq = seq_upper.replace('U', 'T')
        else:
            search_seq = seq_upper.replace('T', 'U')

        matches = re.finditer(pattern, search_seq)

        results = []
        for match in matches:
            start_pos = match.start() + 1
            results.append({
                "start": start_pos,
                "end": match.end(),
                "sequence": match.group()
            })

        return {
            "motif_pattern": pattern,
            "total_matches": len(results),
            "matches": results
        }

    except re.error:
        return {"ERROR": "Invalid motif pattern (regular expression error)."}
    except Exception:
        return {"ERROR": "Unexpected error during motif search."}

# ----------------- 11. K-mer Analysis -----------------
def kmer_analysis(sequence: str, k_length: int) -> dict:
    """
    Performs K-mer counting and returns the frequency of each unique k-mer.
    """
    if k_length <= 0 or k_length > len(sequence):
        return {"ERROR": "K-mer length is invalid or longer than the sequence."}

    kmer_counts = {}
    seq_upper = sequence.upper()
    n = len(seq_upper)

    for i in range(n - k_length + 1):
        kmer = seq_upper[i:i + k_length]
        if 'N' in kmer:
            continue
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    total_kmers = sum(kmer_counts.values())

    kmer_frequencies = {}
    for kmer, count in kmer_counts.items():
        frequency = round((count / total_kmers) * 100, 4) if total_kmers > 0 else 0.0
        kmer_frequencies[kmer] = {
            "count": count,
            "frequency_percent": frequency
        }

    return {
        "k_length": k_length,
        "total_unique_kmers": len(kmer_counts),
        "total_kmers_counted": total_kmers,
        "kmer_data": kmer_frequencies
    }

# ----------------- 12. Sequence Similarity Search (BLAST-like) -----------------
def sequence_similarity_search(query_seq: str, target_seq: str, alignment_type: str = 'local') -> dict:
    """
    Performs a local or global sequence alignment between two sequences
    using Bio.pairwise2 (Smith-Waterman or Needleman-Wunsch).
    """
    if not query_seq or not target_seq:
        return {"ERROR": "Query and target sequences cannot be empty."}

    query_upper = query_seq.upper()
    target_upper = target_seq.upper()

    match_score = 2
    mismatch_score = -1
    gap_open_penalty = -0.5
    gap_extend_penalty = -0.1

    try:
        if alignment_type == 'local':
            alignments = pairwise2.align.localms(
                query_upper, target_upper, match_score, mismatch_score,
                gap_open_penalty, gap_extend_penalty
            )
        else:
            alignments = pairwise2.align.globalms(
                query_upper, target_upper, match_score, mismatch_score,
                gap_open_penalty, gap_extend_penalty
            )

        if not alignments:
            return {"ERROR": "No significant alignment found."}

        best_alignment = alignments[0]
        aligned_query, aligned_target, score, start, end = best_alignment

        identity = sum(1 for a, b in zip(aligned_query, aligned_target) if a == b and a != '-')
        total_aligned_length = len(aligned_query)
        percent_identity = round((identity / total_aligned_length) * 100, 2)

        return {
            "score": round(score, 2),
            "percent_identity": percent_identity,
            "aligned_query": aligned_query,
            "aligned_target": aligned_target,
            "alignment_type": alignment_type,
            "start_target_index": start + 1,
            "end_target_index": end
        }

    except Exception as e:
        return {"ERROR": f"Alignment calculation failed: {str(e)}"}

# ----------------- 13. RNA Secondary Structure Prediction -----------------
def predict_rna_structure(rna_sequence: str) -> dict:
    """
    Predicts structural stability by calculating Tm for a theoretical helix,
    and returns a simple structure representation (Simulated).
    """
    rna_upper = rna_sequence.upper().replace('T', 'U')
    if len(rna_upper) < 15:
        return {"ERROR": "Sequence must be at least 15 bases for stability analysis."}

    helix_seq = rna_upper[:10]

    try:
        tm_stability = round(mt.Tm_GC(helix_seq, seq_type='RNA'), 2)
    except Exception:
        tm_stability = 0.0

    structure_mock = "(((((....)))))"

    return {
        "stability_tm": tm_stability,
        "structure_dot_bracket": structure_mock,
        "sequence_length": len(rna_upper),
        "note": "Structure is a placeholder; Tm is calculated stability."
    }

# ----------------- 14. Signal Peptide Predictor -----------------
def predict_signal_peptide(protein_sequence: str) -> dict:
    """
    Predicts signal peptide region and cleavage site based on basic hydrophobicity rules (Simulated).
    Requires a protein sequence (one-letter amino acid code).
    """
    protein_upper = protein_sequence.upper().replace('*', '')
    if len(protein_upper) < 20:
        return {"ERROR": "Protein sequence must be at least 20 amino acids long."}

    try:
        hydrophobic_core_pattern = r'[LVIAFPMG]{8,16}'
        match = re.search(hydrophobic_core_pattern, protein_upper)

        if match:
            h_region_end_index = match.end()
            cleavage_site = h_region_end_index + 2
            if 15 < cleavage_site < 35:
                potential_signal_peptide = protein_upper[:cleavage_site]
                analyzer = ProteinAnalysis(potential_signal_peptide)
                return {
                    "prediction": "Signal Peptide Detected",
                    "cleavage_site_index": cleavage_site,
                    "signal_peptide_sequence": potential_signal_peptide,
                    "is_hydrophobic": True,
                    "avg_hydrophobicity": round(analyzer.gravy(), 3),
                    "note": "Prediction based on basic hydrophobic core detection rules."
                }

        return {"prediction": "No Signal Peptide Detected", "note": "Failed basic hydrophobic core test."}

    except Exception as e:
        return {"ERROR": f"Prediction failed: {str(e)}"}

# ----------------- 15. Transmembrane Domain Finder -----------------
def find_transmembrane_domains(protein_sequence: str) -> dict:
    """
    Finds transmembrane domains based on a simple sliding window hydrophobicity threshold.
    """
    protein_upper = protein_sequence.upper().replace('*', '')
    if len(protein_upper) < 30:
        return {"ERROR": "Protein sequence must be at least 30 amino acids long."}

    TM_WINDOW = 21
    HYDROPHOBICITY_THRESHOLD = 1.5

    domains = []

    for i in range(len(protein_upper) - TM_WINDOW + 1):
        segment = protein_upper[i:i + TM_WINDOW]
        try:
            analyzer = ProteinAnalysis(segment)
            gravy_score = analyzer.gravy()
            if gravy_score > HYDROPHOBICITY_THRESHOLD:
                domains.append({
                    "start": i + 1,
                    "end": i + TM_WINDOW,
                    "sequence": segment,
                    "hydrophobicity_gravy": round(gravy_score, 3)
                })
        except Exception:
            continue

    if not domains:
        return {"ERROR": "No transmembrane domains found above the GRAVY threshold of 1.5."}

    return {
        "total_domains_found": len(domains),
        "domains": domains
    }
