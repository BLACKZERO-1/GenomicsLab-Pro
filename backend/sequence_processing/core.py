from Bio.SeqUtils import gc_fraction 
from Bio.Seq import Seq
# Import Bio.Restriction module
from Bio.Restriction import RestrictionBatch, EcoRI, BamHI, HindIII, NotI
from Bio.SeqUtils import MeltingTemp as mt 
import numpy as np 

# --- REMOVED: get_codon_usage_data helper function definition from here ---

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
    Uses Bio.Seq.translate().
    """
    try:
        rna_seq_obj = Seq(rna_sequence.upper())
        protein_seq_obj = rna_seq_obj.translate(
            table=1,             # Standard genetic code
            stop_symbol="*",     # Use '*' for stop codons
            to_stop=False        # Translate beyond the first stop codon
        )
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
        # This function calls calculate_gc_content (internal dependency)
        gc_perc = calculate_gc_content(str(window_seq)) 
        gc_list.append(gc_perc)
        
    return gc_list

# ----------------- 7. Codon Usage Analysis -----------------
def analyze_codon_usage(rna_sequence: str) -> dict:
    """
    Calculates the usage frequency (per 1000 codons) and total count 
    for each codon in a coding sequence.
    """
    # Helper function is defined LOCALLY to prevent circular import issues
    def get_codon_usage_data(seq: str) -> dict:
        codon_counts = {}
        for i in range(0, len(seq) - len(seq) % 3, 3):
            codon = seq[i:i+3]
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
        return codon_counts

    try:
        rna_seq_upper = rna_sequence.upper()
        if len(rna_seq_upper) % 3 != 0:
            raise ValueError("Sequence length must be a multiple of 3 for codon analysis.")
        
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

    except ValueError as e:
        return {"ERROR": str(e)}
    except Exception:
        return {"ERROR": "Invalid RNA sequence for codon usage analysis."}

# ----------------- 8. Restriction Site Finder -----------------
def find_restriction_sites(dna_sequence: str) -> dict:
    """
    Finds cutting sites for a standard batch of restriction enzymes 
    (e.g., EcoRI, BamHI, HindIII) within the DNA sequence.
    """
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
        
        return result
        
    except ValueError as e:
        return {"ERROR": f"Validation Error: {str(e)}"}
    except Exception as e:
        print(f"Restriction Site Analysis Failed: {e}") 
        return {"ERROR": "Unexpected error during restriction analysis. Check sequence format."}

# ----------------- 9. Primer Design Tool -----------------
def design_pcr_primers(target_sequence: str, desired_tm: float = 60.0) -> dict:
    """
    Performs basic validation and estimates properties for a theoretical primer pair.
    (A full primer design tool requires complex optimization which is simplified here).
    """
    if len(target_sequence) < 50:
        return {"ERROR": "Target sequence must be at least 50 bases for meaningful primer design."}

    # Use a fixed short sequence from the start for Forward Primer (18 bases)
    forward_seq = target_sequence[:18].upper()
    
    # Use a fixed short sequence from the end for Reverse Primer (18 bases)
    # The reverse primer is the reverse complement of the target's 3' end
    reverse_target = target_sequence[-18:].upper()
    reverse_seq = str(Seq(reverse_target).reverse_complement())
    
    # Calculate properties
    try:
        # Use the imported mt (MeltingTemp) alias from Bio.SeqUtils
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
