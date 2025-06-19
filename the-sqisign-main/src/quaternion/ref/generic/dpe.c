#include "../internal_quaternion_headers/dpe.h"
#include "intbig.h"
#include <math.h>

void
dpe_set_z(dpe_t x, const ibz_t *y) // 포인터로 변경
{
    if (ibz_is_zero(y)) {
        DPE_MANT(x) = 0.0;
        DPE_EXP(x) = DPE_EXPMIN;
        return;
    }

    int is_neg = ibz_is_negative(y);
    ibz_t abs_y;
    if (is_neg) {
        ibz_neg(&abs_y, y);
    } else {
        ibz_copy(&abs_y, y);
    }

    int bits = ibz_bitsize(&abs_y);

    if (bits > 53) {
        ibz_t shifted;
        ibz_div_2exp(&shifted, &abs_y, bits - 53);
        DPE_MANT(x) = (double)shifted[0];
        DPE_EXP(x) = bits - 53;
    } else {
        DPE_MANT(x) = (double)abs_y[0];
        DPE_EXP(x) = 0;
    }

    if (is_neg) {
        DPE_MANT(x) = -DPE_MANT(x);
    }

    dpe_normalize(x);
}

void
dpe_get_z(ibz_t *x, const dpe_t y) // 포인터로 변경
{
    long ey = DPE_EXP(y);

    if (DPE_MANT(y) == 0.0) {
        ibz_init(x);
        return;
    }

    if (ey >= DPE_BITSIZE) {
        double d = DPE_MANT(y) * DPE_2_POW_BITSIZE;
        ibz_init(x);
        (*x)[0] = (uint64_t)fabs(d);
        if (ey > DPE_BITSIZE) {
            ibz_mul_2exp(x, x, (unsigned long)(ey - DPE_BITSIZE));
        }
    } else if (ey < 0) {
        ibz_init(x);
        return;
    } else {
        double d = ldexp(DPE_MANT(y), ey);
        ibz_init(x);
        (*x)[0] = (uint64_t)fabs(round(d));
    }

    if (DPE_MANT(y) < 0.0 && !ibz_is_zero(x)) {
        ibz_neg(x, x);
    }
}