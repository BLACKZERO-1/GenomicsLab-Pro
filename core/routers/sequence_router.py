from fastapi import APIRouter
from pydantic import BaseModel, Field
# Import ALL nine logic functions from the backend file
from backend.sequence_processing.core import (
    calculate_gc_content, 
    transcribe_dna_to_rna, 
    translate_rna_to_protein, 
    get_reverse_complement, 
    find_orfs,
    gc_sliding_window,
    analyze_codon_usage,
    find_restriction_sites,
    design_pcr_primers  # <-- NEW IMPORT
) 

router = APIRouter()

# --- SHARED REQUEST MODEL ---
class SequenceInput(BaseModel):
    sequence: str

# -------------------------- 1. GC Content Endpoint --------------------------
class GCContentResponse(BaseModel):
    sequence_length: int
    gc_percentage: float

@router.post("/gc_content", response_model=GCContentResponse, tags=["Sequence Analysis"])
async def run_gc_content(request: SequenceInput):
    """API endpoint for the GC Content Calculator tool."""
    seq = request.sequence.strip()
    gc_perc = calculate_gc_content(seq)
    return GCContentResponse(sequence_length=len(seq), gc_percentage=gc_perc)

# -------------------------- 2. Transcription Endpoint --------------------------
class TranscriptionResponse(BaseModel):
    input_dna: str
    output_rna: str
    length: int

@router.post("/transcribe", response_model=TranscriptionResponse, tags=["Sequence Analysis"])
async def run_transcription(request: SequenceInput):
    """API endpoint for DNA to RNA Transcription (T -> U)."""
    dna_seq = request.sequence.strip()
    rna_seq = transcribe_dna_to_rna(dna_seq)
    return TranscriptionResponse(input_dna=dna_seq, output_rna=rna_seq, length=len(rna_seq))

# -------------------------- 3. Translation Endpoint --------------------------
class TranslationResponse(BaseModel):
    input_rna: str
    output_protein: str
    length: int

@router.post("/translate", response_model=TranslationResponse, tags=["Sequence Analysis"])
async def run_translation(request: SequenceInput):
    """API endpoint for RNA to Protein Translation (codon to amino acid)."""
    rna_seq = request.sequence.strip()
    protein_seq = translate_rna_to_protein(rna_seq)
    return TranslationResponse(input_rna=rna_seq, output_protein=protein_seq, length=len(protein_seq))

# -------------------------- 4. Reverse Complement Endpoint --------------------------
class ReverseComplementResponse(BaseModel):
    input_sequence: str
    reverse_complement: str
    length: int

@router.post("/reverse_complement", response_model=ReverseComplementResponse, tags=["Sequence Analysis"])
async def run_reverse_complement(request: SequenceInput):
    """API endpoint to calculate the reverse complement of a DNA sequence."""
    dna_seq = request.sequence.strip()
    rev_comp = get_reverse_complement(dna_seq)
    return ReverseComplementResponse(input_sequence=dna_seq, reverse_complement=rev_comp, length=len(rev_comp))

# -------------------------- 5. ORF Finder Endpoint --------------------------
class ORFResponse(BaseModel):
    orfs_found: dict = Field(
        example={"Frame +1": ["MSRF..."]},
        description="Dictionary of identified peptides (strings) grouped by reading frame."
    )
    total_frames_searched: int = 6

@router.post("/orf_finder", response_model=ORFResponse, tags=["Sequence Analysis"])
async def run_orf_finder(request: SequenceInput):
    """API endpoint to find Open Reading Frames (ORFs) in all six frames."""
    dna_seq = request.sequence.strip()
    orfs = find_orfs(dna_seq)
    return ORFResponse(orfs_found=orfs, total_frames_searched=6)


# -------------------------- 6. GC Sliding Window Endpoint --------------------------
class GCSlidingWindowRequest(BaseModel):
    sequence: str = Field(..., description="The input DNA sequence.")
    window_size: int = Field(100, ge=10, description="The size of the sliding window.")
    step_size: int = Field(50, ge=1, description="The distance to slide the window.")

class GCSlidingWindowResponse(BaseModel):
    window_size: int
    step_size: int
    gc_percentages: list[float] = Field(..., description="List of GC percentages for each window.")

@router.post("/gc_sliding_window", response_model=GCSlidingWindowResponse, tags=["Sequence Analysis"])
async def run_gc_sliding_window(request: GCSlidingWindowRequest):
    """API endpoint to calculate GC content across a sequence using a sliding window."""
    seq = request.sequence.strip()
    gc_list = gc_sliding_window(seq, request.window_size, request.step_size)
    
    return GCSlidingWindowResponse(
        window_size=request.window_size,
        step_size=request.step_size,
        gc_percentages=gc_list
    )


# -------------------------- 7. Codon Usage Analysis Endpoint --------------------------
class CodonUsageResponse(BaseModel):
    total_codons: int
    codon_usage: dict = Field(
        example={
            "AUG": {"count": 10, "frequency_per_1000": 30.3},
            "UAC": {"count": 15, "frequency_per_1000": 45.45}
        },
        description="Usage count and frequency per 1000 for each codon found."
    )

@router.post("/codon_usage", response_model=CodonUsageResponse, tags=["Sequence Analysis"])
async def run_codon_usage(request: SequenceInput):
    """API endpoint for Codon Usage Analysis of a coding RNA sequence."""
    rna_seq = request.sequence.strip()
    usage_data = analyze_codon_usage(rna_seq)
    
    if "ERROR" in usage_data:
        return CodonUsageResponse(total_codons=0, codon_usage=usage_data)

    return CodonUsageResponse(
        total_codons=usage_data['total_codons'],
        codon_usage=usage_data['codon_usage']
    )

# -------------------------- 8. Restriction Site Finder Endpoint --------------------------
class RestrictionSiteResponse(BaseModel):
    enzymes_tested: list[str] = ["EcoRI", "BamHI", "HindIII", "NotI"]
    sites_found: dict = Field(
        example={
            "EcoRI": [105, 520],
            "BamHI": [350]
        },
        description="Dictionary mapping enzyme name to a list of 1-based cutting positions."
    )

@router.post("/restriction_sites", response_model=RestrictionSiteResponse, tags=["Sequence Analysis"])
async def run_restriction_finder(request: SequenceInput):
    """API endpoint to find cutting sites for common restriction enzymes."""
    dna_seq = request.sequence.strip()
    sites = find_restriction_sites(dna_seq)
    
    if "ERROR" in sites:
        return RestrictionSiteResponse(sites_found={"ERROR": sites["ERROR"]}, enzymes_tested=[])

    tested_enzymes = list(sites.keys())
    
    return RestrictionSiteResponse(
        enzymes_tested=tested_enzymes,
        sites_found=sites
    )

# -------------------------- 9. Primer Design Endpoint --------------------------
class PrimerDesignRequest(BaseModel):
    target_sequence: str = Field(..., min_length=50, description="The DNA sequence to amplify (must be >= 50 bases).")
    desired_tm: float = Field(60.0, ge=50.0, le=75.0, description="The desired target melting temperature in Celsius.")

class PrimerDesignResponse(BaseModel):
    Summary: dict
    Forward_Primer: dict
    Reverse_Primer: dict

@router.post("/primer_design", response_model=PrimerDesignResponse, tags=["Sequence Analysis"])
async def run_primer_design(request: PrimerDesignRequest):
    """API endpoint for basic PCR Primer Design (fixed 18-mer primers)."""
    # Note: The request uses target_sequence and desired_tm, but the logic only uses target_sequence
    result = design_pcr_primers(request.target_sequence.strip(), request.desired_tm)
    
    if "ERROR" in result:
        return PrimerDesignResponse(
            Summary={"ERROR": result["ERROR"]},
            Forward_Primer={},
            Reverse_Primer={}
        )
        
    return PrimerDesignResponse(
        Summary=result["Summary"],
        Forward_Primer=result["Forward_Primer"],
        Reverse_Primer=result["Reverse_Primer"]
    )
