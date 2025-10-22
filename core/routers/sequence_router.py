# core/routers/sequence_router.py
import asyncio
from fastapi import APIRouter
from pydantic import BaseModel, Field
from typing import List, Dict

# Import logic functions from backend core
from backend.sequence_processing.core import (
    calculate_gc_content,
    transcribe_dna_to_rna,
    translate_rna_to_protein,
    get_reverse_complement,
    find_orfs,
    gc_sliding_window,
    analyze_codon_usage,
    find_restriction_sites,
    design_pcr_primers,
    find_sequence_motifs,
    kmer_analysis,
    sequence_similarity_search,
    predict_rna_structure,
    predict_signal_peptide,
    find_transmembrane_domains
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
    gc_perc = await asyncio.to_thread(calculate_gc_content, seq)
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
    rna_seq = await asyncio.to_thread(transcribe_dna_to_rna, dna_seq)
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
    protein_seq = await asyncio.to_thread(translate_rna_to_protein, rna_seq)
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
    rev_comp = await asyncio.to_thread(get_reverse_complement, dna_seq)
    return ReverseComplementResponse(input_sequence=dna_seq, reverse_complement=rev_comp, length=len(rev_comp))

# -------------------------- 5. ORF Finder Endpoint --------------------------
class ORFResponse(BaseModel):
    orfs_found: Dict[str, List[str]] = Field(example={"Frame +1": ["MSRF..."]}, description="Dictionary of identified peptides grouped by reading frame.")
    total_frames_searched: int = 6

@router.post("/orf_finder", response_model=ORFResponse, tags=["Sequence Analysis"])
async def run_orf_finder(request: SequenceInput):
    """API endpoint to find Open Reading Frames (ORFs) in all six frames."""
    dna_seq = request.sequence.strip()
    orfs = await asyncio.to_thread(find_orfs, dna_seq)
    return ORFResponse(orfs_found=orfs, total_frames_searched=6)


# -------------------------- 6. GC Sliding Window Endpoint --------------------------
class GCSlidingWindowRequest(BaseModel):
    sequence: str = Field(..., description="The input DNA sequence.")
    window_size: int = Field(100, ge=10, description="The size of the sliding window.")
    step_size: int = Field(50, ge=1, description="The distance to slide the window.")

class GCSlidingWindowResponse(BaseModel):
    window_size: int
    step_size: int
    gc_percentages: List[float] = Field(..., description="List of GC percentages for each window.")

@router.post("/gc_sliding_window", response_model=GCSlidingWindowResponse, tags=["Sequence Analysis"])
async def run_gc_sliding_window(request: GCSlidingWindowRequest):
    """API endpoint to calculate GC content across a sequence using a sliding window."""
    seq = request.sequence.strip()
    gc_list = await asyncio.to_thread(gc_sliding_window, seq, request.window_size, request.step_size)
    return GCSlidingWindowResponse(window_size=request.window_size, step_size=request.step_size, gc_percentages=gc_list)


# -------------------------- 7. Codon Usage Analysis Endpoint --------------------------
class CodonUsageResponse(BaseModel):
    total_codons: int
    codon_usage: Dict[str, dict] = Field(example={"AUG": {"count": 10, "frequency_per_1000": 30.3}}, description="Usage count and frequency per 1000 for each codon found.")

@router.post("/codon_usage", response_model=CodonUsageResponse, tags=["Sequence Analysis"])
async def run_codon_usage(request: SequenceInput):
    """API endpoint for Codon Usage Analysis of a coding RNA sequence."""
    rna_seq = request.sequence.strip()
    usage_data = await asyncio.to_thread(analyze_codon_usage, rna_seq)
    if isinstance(usage_data, dict) and "ERROR" in usage_data:
        return CodonUsageResponse(total_codons=0, codon_usage={})
    return CodonUsageResponse(total_codons=usage_data['total_codons'], codon_usage=usage_data['codon_usage'])

# -------------------------- 8. Restriction Site Finder Endpoint --------------------------
class RestrictionSiteResponse(BaseModel):
    enzymes_tested: List[str] = Field(default_factory=lambda: ["EcoRI", "BamHI", "HindIII", "NotI"])
    sites_found: Dict[str, List[int]] = Field(..., description="Dictionary mapping enzyme name to a list of 1-based cutting positions.")

@router.post("/restriction_sites", response_model=RestrictionSiteResponse, tags=["Sequence Analysis"])
async def run_restriction_finder(request: SequenceInput):
    """API endpoint to find cutting sites for common restriction enzymes."""
    dna_seq = request.sequence.strip()
    sites = await asyncio.to_thread(find_restriction_sites, dna_seq)
    if isinstance(sites, dict) and "ERROR" in sites:
        return RestrictionSiteResponse(enzymes_tested=[], sites_found={})
    tested_enzymes = list(sites.keys())
    return RestrictionSiteResponse(enzymes_tested=tested_enzymes, sites_found=sites)

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
    result = await asyncio.to_thread(design_pcr_primers, request.target_sequence.strip(), request.desired_tm)
    if isinstance(result, dict) and "ERROR" in result:
        return PrimerDesignResponse(Summary={"ERROR": result["ERROR"]}, Forward_Primer={}, Reverse_Primer={})
    return PrimerDesignResponse(Summary=result["Summary"], Forward_Primer=result["Forward_Primer"], Reverse_Primer=result["Reverse_Primer"])

# -------------------------- 10. Motif Search Endpoint --------------------------
class MotifSearchRequest(BaseModel):
    sequence: str = Field(..., description="The sequence to search within (DNA or RNA).")
    motif_pattern: str = Field(..., description="The pattern to search for (e.g., 'ATGG' or ambiguous codes like 'WSRT').")
    is_dna: bool = Field(True, description="Set to True if the sequence is DNA (default).")

class MotifMatch(BaseModel):
    start: int
    end: int
    sequence: str

class MotifSearchResponse(BaseModel):
    motif_pattern: str
    total_matches: int
    matches: List[MotifMatch]

@router.post("/motif_search", response_model=MotifSearchResponse, tags=["Sequence Analysis"])
async def run_motif_search(request: MotifSearchRequest):
    """API endpoint to find sequence motifs using regex/IUPAC codes."""
    result = await asyncio.to_thread(find_sequence_motifs, request.sequence, request.motif_pattern, request.is_dna)
    if isinstance(result, dict) and "ERROR" in result:
        return MotifSearchResponse(motif_pattern=request.motif_pattern, total_matches=0, matches=[])
    return MotifSearchResponse(
        motif_pattern=result["motif_pattern"],
        total_matches=result["total_matches"],
        matches=result["matches"]
    )

# -------------------------- 11. K-mer Analysis Endpoint --------------------------
class KmerAnalysisRequest(BaseModel):
    sequence: str = Field(..., description="The DNA or RNA sequence to analyze.")
    k_length: int = Field(3, ge=2, le=10, description="The length (k) of the oligomer to count (e.g., 3 for trinucleotides).")

class KmerData(BaseModel):
    count: int
    frequency_percent: float

class KmerAnalysisResponse(BaseModel):
    k_length: int
    total_unique_kmers: int
    total_kmers_counted: int
    kmer_data: Dict[str, KmerData] = Field(..., description="Dictionary mapping k-mer sequence to its count and frequency.")

@router.post("/kmer_analysis", response_model=KmerAnalysisResponse, tags=["Sequence Analysis"])
async def run_kmer_analysis(request: KmerAnalysisRequest):
    """API endpoint to perform K-mer frequency analysis."""
    result = await asyncio.to_thread(kmer_analysis, request.sequence.strip(), request.k_length)
    if isinstance(result, dict) and "ERROR" in result:
        return KmerAnalysisResponse(k_length=request.k_length, total_unique_kmers=0, total_kmers_counted=0, kmer_data={})
    # Convert inner structures to Pydantic-friendly form if needed
    return KmerAnalysisResponse(
        k_length=result["k_length"],
        total_unique_kmers=result["total_unique_kmers"],
        total_kmers_counted=result["total_kmers_counted"],
        kmer_data=result["kmer_data"]
    )

# -------------------------- 12. Similarity Search Endpoint --------------------------
class SimilaritySearchRequest(BaseModel):
    query_sequence: str = Field(..., description="The short sequence to search.")
    target_sequence: str = Field(..., description="The long reference sequence to search within.")
    alignment_type: str = Field('local', description="Type of alignment: 'local' (Smith-Waterman) or 'global' (Needleman-Wunsch).")

class SimilaritySearchResponse(BaseModel):
    aligned_query: str
    aligned_target: str
    score: float
    percent_identity: float
    alignment_type: str
    start_target_index: int
    end_target_index: int

@router.post("/similarity_search", response_model=SimilaritySearchResponse, tags=["Sequence Analysis"])
async def run_similarity_search(request: SimilaritySearchRequest):
    """API endpoint for local (BLAST-like) or global sequence comparison."""
    result = await asyncio.to_thread(sequence_similarity_search,
                                     request.query_sequence.strip(),
                                     request.target_sequence.strip(),
                                     request.alignment_type)
    if isinstance(result, dict) and "ERROR" in result:
        return SimilaritySearchResponse(
            aligned_query=result["ERROR"],
            aligned_target="---",
            score=0.0,
            percent_identity=0.0,
            alignment_type=request.alignment_type,
            start_target_index=0,
            end_target_index=0
        )
    return SimilaritySearchResponse(**result)

# -------------------------- 13. RNA Secondary Structure Prediction --------------------------
class SecondaryStructureResponse(BaseModel):
    stability_tm: float = Field(..., description="The calculated stability (Melting Temperature) of the theoretical helix (in Celsius).")
    structure_dot_bracket: str = Field(..., description="Mock secondary structure representation (dot-bracket notation).")
    sequence_length: int
    note: str

@router.post("/rna_structure_prediction", response_model=SecondaryStructureResponse, tags=["Sequence Analysis"])
async def run_rna_structure_prediction(request: SequenceInput):
    """API endpoint for RNA Secondary Structure Prediction (Simulated)."""
    result = await asyncio.to_thread(predict_rna_structure, request.sequence.strip())
    if isinstance(result, dict) and "ERROR" in result:
        return SecondaryStructureResponse(
            stability_tm=0.0,
            structure_dot_bracket="ERROR",
            sequence_length=len(request.sequence.strip()),
            note=result["ERROR"]
        )
    return SecondaryStructureResponse(**result)


# -------------------------- 14. Signal Peptide Predictor Endpoint --------------------------
class SignalPeptideResponse(BaseModel):
    prediction: str = Field(..., description="Prediction status: 'Signal Peptide Detected' or 'No Signal Peptide Detected'.")
    cleavage_site_index: int = Field(0, description="1-based index of the predicted cleavage site, or 0 if none found.")
    signal_peptide_sequence: str = Field("", description="The N-terminal sequence identified as the signal peptide.")
    is_hydrophobic: bool
    avg_hydrophobicity: float
    note: str

@router.post("/signal_peptide_predictor", response_model=SignalPeptideResponse, tags=["Sequence Analysis"])
async def run_signal_peptide_predictor(request: SequenceInput):
    """API endpoint for Signal Peptide and Cleavage Site Prediction (Rule-based Simulation)."""
    result = await asyncio.to_thread(predict_signal_peptide, request.sequence.strip())
    if isinstance(result, dict) and "ERROR" in result:
        return SignalPeptideResponse(
            prediction="ERROR",
            cleavage_site_index=0,
            signal_peptide_sequence=result.get("ERROR", ""),
            is_hydrophobic=False,
            avg_hydrophobicity=0.0,
            note="Prediction failed."
        )
    if result.get("prediction") == "No Signal Peptide Detected":
        return SignalPeptideResponse(
            prediction=result["prediction"],
            cleavage_site_index=0,
            signal_peptide_sequence="",
            is_hydrophobic=False,
            avg_hydrophobicity=0.0,
            note=result.get("note", "")
        )
    return SignalPeptideResponse(
        prediction=result["prediction"],
        cleavage_site_index=result["cleavage_site_index"],
        signal_peptide_sequence=result["signal_peptide_sequence"],
        is_hydrophobic=result["is_hydrophobic"],
        avg_hydrophobicity=result["avg_hydrophobicity"],
        note=result["note"]
    )

# -------------------------- 15. Transmembrane Domain Finder Endpoint --------------------------
class TransmembraneDomain(BaseModel):
    start: int
    end: int
    sequence: str
    hydrophobicity_gravy: float

class TransmembraneDomainResponse(BaseModel):
    total_domains_found: int
    domains: List[TransmembraneDomain] = Field(..., description="List of predicted transmembrane domains and their properties.")

@router.post("/transmembrane_domain_finder", response_model=TransmembraneDomainResponse, tags=["Sequence Analysis"])
async def run_transmembrane_domain_finder(request: SequenceInput):
    """API endpoint to find Transmembrane Domains (Simulated based on hydrophobicity)."""
    result = await asyncio.to_thread(find_transmembrane_domains, request.sequence.strip())
    if isinstance(result, dict) and "ERROR" in result:
        return TransmembraneDomainResponse(
            total_domains_found=0,
            domains=[]
        )
    return TransmembraneDomainResponse(
        total_domains_found=result["total_domains_found"],
        domains=result["domains"]
    )
