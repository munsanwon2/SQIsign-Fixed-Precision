#!/usr/bin/env python3
"""Run the Sage parameter audits and write docs/new-parameter-manifest.json."""

import json
import os
from pathlib import Path
import subprocess

scripts = Path(__file__).resolve().parent
root = scripts.parents[1]
env = os.environ.copy()
env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(scripts), env.get("PYTHONPATH"))))
env["PYTHONDONTWRITEBYTECODE"] = "1"

manifest = {
    "description": "Round-2 algorithms with Round-3 prime/security parameters",
    "baseline": "d0cb037ee6a9f68f55cb6f55b4e3746c79550330",
    "parameter_spec": "https://sqisign.org/spec/sqisign-20260901.pdf",
    "upstream_tag": "nist-v3",
    "upstream_commit": "6d017708db403bf83977fa70770fc4f7f9e9ff21",
    "section3_paper_sha256_reported_in_task": "94d946d798f35798b84aa5de6b69adf2024e0487720f292523d0f6ee250c540c",
    "primality_check": "SageMath is_prime(proof=True) and next_prime(proof=True)",
    "sage": subprocess.check_output(["sage", "--version"], text=True, env=env).strip(),
    "generator_dependencies": {
        "deuring-2D": "799d69d0c784a57446566e15be4207b936bcc469",
        "two-isogenies": "0ea5020effabde4bdc5e6f425ba5b94209f0e58b",
    },
    "levels": [],
}
for level in (1, 3, 5):
    output = subprocess.check_output(
        ["sage", str(scripts / "audit_parameters.sage")],
        text=True, cwd=root / f"src/precomp/ref/lvl{level}", env=env,
    )
    manifest["levels"].append(json.loads(output))
(root / "docs/new-parameter-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
