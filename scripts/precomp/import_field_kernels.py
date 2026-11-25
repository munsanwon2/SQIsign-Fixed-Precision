#!/usr/bin/env python3
"""Import the field kernels from the pinned Apache-2.0 SQIsign submission."""

import argparse
from pathlib import Path
import re
import subprocess

UPSTREAM = "6d017708db403bf83977fa70770fc4f7f9e9ff21"
PRIMES = {1: "p324_3", 3: "p500_27", 5: "p664_17"}
MARKER = "/******************************************************************************\n * SQIsign fp_* API wrappers"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkout", type=Path, help="the-sqisign checkout at nist-v3")
    upstream = parser.parse_args().checkout
    root = Path(__file__).resolve().parents[2] / "src/gf"

    def source(path):
        return subprocess.check_output(
            ["git", "-C", str(upstream), "show", f"{UPSTREAM}:src/gf/{path}"],
            text=True,
        )

    def write(path, content):
        (root / path).write_text(content, encoding="utf-8")

    notice = (
        "// SPDX-License-Identifier: Apache-2.0\n"
        f"// Field kernel from SQIsign nist-v3 ({UPSTREAM}).\n"
    )
    for level, prime in PRIMES.items():
        for word_size in (32, 64):
            filename = f"fp_{prime}_{word_size}.inc"
            body, suffix = source(f"ref/{prime}/{filename}").split(MARKER)
            assert "fp_decode(" in suffix
            write(f"ref/lvl{level}/{filename}", notice + body.rstrip() + "\n")

        asm = source(f"broadwell/{prime}/fp_asm.S")
        # Bind the external modulus symbols through the local namespace.
        asm = re.sub(r"\[rip\+(p2?)\b",
                     lambda m: "[rip+" + {"p": "fp_modulus", "p2": "fp_twice_modulus"}[m[1]], asm)
        write(f"broadwell/lvl{level}/fp_asm.S", notice + asm)


if __name__ == "__main__":
    main()
