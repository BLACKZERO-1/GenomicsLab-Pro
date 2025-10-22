# core/app.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from core.routers.sequence_router import router as sequence_router

app = FastAPI(title="GenomicsLab Pro API", version="1.0.0")

# --- CORS MIDDLEWARE (development) ---
origins = ["*"]  # For development allow all; lock down in production
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
# --------------------------------------

app.include_router(sequence_router, prefix="/sequence", tags=["Sequence Analysis"])

@app.get("/")
def read_root():
    """Root endpoint to confirm the API is running."""
    return {"Welcome": "GenomicsLab Pro - Backend Active"}

@app.get("/health")
def health_check():
    """Simple endpoint for frontend connection check."""
    return {"status": "ok"}

@app.on_event("startup")
def on_startup():
    print("🚀 GenomicsLab Pro backend started successfully.")
