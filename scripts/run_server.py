# runserver.py
import uvicorn
import os
import sys

# --- Ensure project root is in path ---
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
# ----------------------------------------------------

if __name__ == "__main__":
    print(f"Starting GenomicsLabPro backend from: {PROJECT_ROOT}")
    uvicorn.run("core.app:app", host="127.0.0.1", port=8000, reload=True)
