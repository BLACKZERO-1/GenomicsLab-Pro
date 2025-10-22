from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware # <--- NEW IMPORT
from core.routers.sequence_router import router as sequence_router 
    
app = FastAPI(title="GenomicsLab Pro API", version="1.0.0")
    
# --- CRITICAL FIX: CORS MIDDLEWARE ---
origins = [
    "*", # Allow all origins for local development simplicity
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"], # Allow all HTTP methods (POST, GET)
    allow_headers=["*"], # Allow all headers
)
# ------------------------------------

app.include_router(sequence_router, prefix="/sequence", tags=["Sequence Analysis"]) 
    
@app.get("/")
def read_root():
    """Root endpoint to confirm the API is running."""
    return {"Welcome": "GenomicsLab Pro - Backend Active"}
