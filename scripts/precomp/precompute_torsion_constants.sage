#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import f, ibz_limbs, mixing_degree, security_bits, torsion_bytes

################################################################

tors2part = 1 << f
deg_mix = mixing_degree()

defs = {
        'TORSION_2POWER_BYTES': torsion_bytes,
    }

from cformat import Ibz, Object, ObjectFormatter

objs = ObjectFormatter([
        Object('ibz_t', 'TWO_TO_SECURITY_BITS', Ibz(1 << security_bits)),
        Object('ibz_t', 'TORSION_PLUS_2POWER', Ibz(tors2part)),
        Object('ibz_t', 'SEC_DEGREE', Ibz(deg_mix)),
        Object('ibz_t', 'COM_DEGREE', Ibz(deg_mix)),
    ])

with open('include/torsion_constants.h','w') as hfile:
    with open('torsion_constants.c','w') as cfile:
        print(f'#include <quaternion.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <torsion_constants.h>', file=cfile)
        print(f'\n_Static_assert(IBZ_LIMBS == {ibz_limbs}, "precomputed ibz_t width mismatch");\n', file=cfile)

        for k,v in defs.items():
            print(f'#define {k} {v}', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

