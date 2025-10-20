import uvicorn
import os
import sys

# --- CRITICAL FIX: Ensure project root is in path ---
# 1. Gets the absolute path to the GenomicsLabPro root folder.
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# 2. Inserts the project root path at the beginning of Python's import search path.
# This guarantees modules like 'core' and 'backend' are found instantly.
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
# ----------------------------------------------------

if __name__ == "__main__":
    # Runs the app instance called 'app' inside the 'core.app' module.
    uvicorn.run("core.app:app", host="127.0.0.1", port=8000, reload=True)
