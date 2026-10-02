"""Set a writable XML cache, then execute the selected refresh script."""

import os
from pathlib import Path
import runpy
import sys

from Bio import Entrez

Entrez.local_cache = os.environ["GPCR_PENDB_REFRESH_XML_CACHE"]
script = sys.argv[1]
sys.argv = sys.argv[1:]
sys.path.insert(0, str(Path(script).resolve().parent))
runpy.run_path(script, run_name="__main__")
