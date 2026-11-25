#!/usr/bin/env sage
"""Check a level's inputs and generated integer constants, then emit its manifest."""

import json
import re
from pathlib import Path

from cformat import Ibz
from parameters import (p, p_bits, f, cofactor, lvl, mixing_degree, security_bits,
                        num_orders, ibz_limbs, ref_radix, fp_bytes, response_length,
                        challenge_space, challenge_bits, hash_iterations, response_bytes,
                        counter_bytes, secret_ideal_bytes, torsion_bytes, hd_extra_torsion)

expected = {
    1: (3, 324, 128, 582, 2043, 9555, 9565, 150),
    3: (27, 500, 192, 889, 3118, 14572, 14569, 228),
    5: (17, 664, 256, 1181, 4098, 19145, 19145, 300),
}
c, exponent, security, mix_bits, b_keygen, b_paper, b_line13, limbs = expected[lvl]
assert p == (c << exponent) - 1 and cofactor == c and f == exponent
assert security_bits == security and p % 4 == 3 and p.is_prime(proof=True)

mix = mixing_degree()  # Sage proves primality of the first prime above the threshold.
assert mix.nbits() == mix_bits
assert ibz_limbs == limbs

RIF = RealIntervalField(256)
def bits(bound):
    log2 = bound.log() / RIF(2).log()
    b = ZZ(log2.lower().ceil())
    assert log2.upper().ceil() == b
    return b

P = RIF(p)
C = 47 / RIF.pi()**2
keygen = bits(2**80 * C**4 * P**6)
paper = bits(2**384 * C**24 * P**28)
line13 = bits(2**16 * (2**16 * C * P)**24 * 2**(4*f))
assert (keygen, paper, line13) == (b_keygen, b_paper, b_line13)
storage_bits = 64 * ibz_limbs
capacity = storage_bits - 1
assert capacity >= max(paper, line13)
assert 64 * (ibz_limbs - 1) - 1 < max(paper, line13)

def decode_ibz(words):
    assert len(words) == ibz_limbs
    value = sum(word << (64 * i) for i, word in enumerate(words))
    if value >> capacity:
        value -= 1 << storage_bits
    return value

limit = 1 << capacity
for value in (0, 1, -1, 2**64, -(2**64 + 1), limit - 1, -limit):
    assert decode_ibz(Ibz(value).limbs()) == value
for value in (limit, -limit - 1):
    try:
        Ibz(value).limbs()
    except OverflowError:
        pass
    else:
        raise AssertionError("out-of-range initializer")

tokens = re.compile(r"UINT64_C\((0x[0-9a-f]+)\)|UINT64_MAX|\b0\b")
def decode_literal(literal):
    words = []
    for token in tokens.finditer(literal):
        if token[1]:
            words.append(int(token[1], 16))
        else:
            words.append(2**64 - 1 if token[0] == "UINT64_MAX" else 0)
    return decode_ibz(words)

constants = Path("torsion_constants.c").read_text()
for symbol, value in {
    "TWO_TO_SECURITY_BITS": 1 << security_bits, "TORSION_PLUS_2POWER": 1 << f,
    "SEC_DEGREE": mix, "COM_DEGREE": mix,
}.items():
    literal = re.search(rf"const ibz_t {symbol} = (\{{[^}}]*\}});", constants)[1]
    assert decode_literal(literal) == value

proof.all(False)
from maxorders import orders, Quat1

order_q = [q for q, *_ in orders]
connecting_norms = [Quat1.ideal(list(map(Quat1, data[4]))).norm() for data in orders]
assert len(order_q) == num_orders and max(order_q) <= 181
assert all(norm % 2 for norm in connecting_norms)
assert all(norm**2 < q * p for norm, q in zip(connecting_norms, order_q))

quat_data = Path("quaternion_data.c").read_text()
negative_literals = [decode_literal(m[0]) for m in re.finditer(r"\{\n    [^{}]*?\n\}", quat_data)]
assert any(value < 0 for value in negative_literals)
for _, _, mat, _, ideal, gamma in orders:
    for data in (mat, ideal, gamma):
        for value in (data * data.denominator()).list():
            if value < 0:
                assert Ibz(value)._literal() in quat_data

def field_layout(radix):
    n = (p_bits + radix - 1) // radix
    r = 1 << (radix * n)
    return dict(radix_bits=radix, limbs=n, montgomery_r_bits=radix*n,
                montgomery_r_mod_p_hex=hex(r % p), montgomery_r2_mod_p_hex=hex(r**2 % p),
                montgomery_reduction_inverse=int((-inverse_mod(p, 1 << radix)) % (1 << radix)))

manifest = dict(
    level=int(lvl), p_hex=hex(p), p_bits=int(p_bits), c=int(cofactor), f=int(f),
    security_bits=int(security_bits), mixing_degree_hex=hex(mix), mixing_degree_bits=int(mix_bits),
    response_exponent=int(response_length), challenge_space_bits=int(challenge_space),
    challenge_bits=int(challenge_bits), hash_iterations=int(hash_iterations),
    hd_extra_torsion=int(hd_extra_torsion),
    paper_keygen_bits=int(keygen), paper_reported_bits=int(paper), line13_expression_bits=int(line13),
    verified_section3_expression_bits=int(max(paper, line13)), validated_required_bits=None,
    ibz_limbs=int(ibz_limbs), ibz_bits=int(storage_bits), signed_capacity_bits=int(capacity),
    field_bytes=int(fp_bytes), scalar_storage_bytes=int(fp_bytes), torsion_integer_bytes=int(torsion_bytes),
    response_matrix_entry_bytes=int(response_bytes), response_counter_bytes=int(counter_bytes),
    secret_ideal_entry_bytes=int(secret_ideal_bytes), public_key_bytes=int(2*fp_bytes + 1),
    secret_key_bytes=int(2*fp_bytes + 1 + 5*secret_ideal_bytes + 4*torsion_bytes),
    signature_bytes=int(2*fp_bytes + 2*counter_bytes + 4*response_bytes + security_bits//8 + 2),
    field_layouts={
        "ref64": field_layout(ref_radix[64]),
        "ref32": field_layout(ref_radix[32]),
        "broadwell": field_layout(64),
    },
    extremal_order_q=list(map(int, order_q)), connecting_ideal_norms_hex=list(map(hex, connecting_norms)),
)
print(json.dumps(manifest, default=int))
