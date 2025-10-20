from fastapi import FastAPI
# Import the sequence router which contains the GC Content, Transcription, and Translation endpoints.
from core.routers.sequence_router import router as sequence_router 
    
app = FastAPI(title="GenomicsLab Pro API", version="1.0.0")
    
# Link the sequence router. All endpoints in it will start with the prefix /sequence.
app.include_router(sequence_router, prefix="/sequence", tags=["Sequence Analysis"]) 
    
@app.get("/")
def read_root():
    """Root endpoint to confirm the API is running."""
    return {"Welcome": "GenomicsLab Pro - Backend Active"}
