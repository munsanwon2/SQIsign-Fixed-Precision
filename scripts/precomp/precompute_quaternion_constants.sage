#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import p_bits, f
negl = 2**-64

################################################################

defs = dict()

# RepresentInteger data
small = ceil(log(negl, 2) / -1)
assert 2**-small <= negl

ri = RealIntervalField(128)
failure_probability = ri(1 - 1 / (64 * p_bits))
shift = (ri(negl).log() / failure_probability.log()).log() / ri(2).log()
add_shift = ZZ(shift.upper().ceil())
assert shift.lower().ceil() == add_shift

defs['QUAT_primality_num_iter'] = ceil(-log(negl, 4))
defs['QUAT_repres_bound_input'] = add_shift

# Equivalent ideal data
defs['QUAT_equiv_bound_coeff'] = 2**(1 + add_shift//4)

# Find_uv constants
m = 2 + (p_bits - f) // 4
defs['FINDUV_box_size'] = m
defs['FINDUV_cube_size'] = (2 * m + 1)**4 - 1

################################################################

with open('include/quaternion_constants.h','w') as hfile:
    print(f'#include <quaternion.h>', file=hfile)

    for k,v in defs.items():
        v = ZZ(v)
        print(f'#define {k} {v}', file=hfile)


