"""Require the frozen runner and block network use in every packaged code file."""

import os
from pathlib import Path
import socket
import sys


ROOT = Path(__file__).resolve().parent.parent
GUARD = ROOT / "Offline_Guard"


def _network_blocked(*args, **kwargs):
    raise SystemExit("FROZEN WORKFLOW: network access is disabled")


socket.socket.connect = _network_blocked
socket.socket.connect_ex = _network_blocked
socket.create_connection = _network_blocked
sys.path.insert(0, str(GUARD))

if os.environ.get("GPCR_PENDB_FROZEN_RUNNER") != str(ROOT):
    raise SystemExit(
        "FROZEN WORKFLOW: run each stage with run_frozen.py --work DIR "
        "--until S<N>. Direct Code S1-S15 invocation is disabled so it "
        "cannot query live sources or read a different local input."
    )
