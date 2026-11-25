#!/usr/bin/env python3
from pathlib import Path

from sage.all import ZZ

inputs = {}
for line in Path("sqisign_parameters.txt").read_text().splitlines():
    line = line.split("#", 1)[0].strip()
    if not line:
        continue
    key, value = (part.strip() for part in line.split("="))
    if key in inputs:
        raise ValueError(f"duplicate parameter {key}")
    inputs[key] = ZZ(value, 0)

lvl = inputs["lvl"]
p = inputs["p"]
security_bits = inputs["security_bits"]
num_orders = inputs["num_orders"]
ibz_limbs = inputs["ibz_limbs"]
ref_radix = {32: inputs["ref_radix32"], 64: inputs["ref_radix64"]}

p_bits = p.nbits()
f = (p + 1).valuation(2)
cofactor = (p + 1) >> f
assert p % 4 == 3 and f > 0
assert security_bits > 0 and security_bits % 8 == 0

response_length = (p_bits + 1) // 2
challenge_space = f - response_length
assert challenge_space > 0
challenge_bits = min(security_bits, challenge_space)
hash_iterations = 1 << (security_bits - challenge_bits)

fp_bytes = ((p_bits + 63) // 64) * 8
torsion_bytes = (f + 8) // 8  # Includes the bit for the value 2^f.
hd_extra_torsion = 2
response_bytes = (response_length + hd_extra_torsion + 7) // 8
counter_bytes = (response_length.bit_length() + 7) // 8
assert counter_bytes in (1, 2)
assert response_length + hd_extra_torsion <= f

# Lemma 16 places the equivalent prime norm below 2^19 * p.
# The generator needs at most two HNF columns (each bounded by 2 * norm).
secret_ideal_bytes = max(fp_bytes, (p_bits + 19 + 2 + 1 + 7) // 8)

def mixing_degree():
    return (p << (2 * security_bits)).next_prime(proof=True)

