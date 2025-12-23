#include "../internal_quaternion_headers/dpe.h"
#include "intbig.h"
#include <math.h>

static void ibz_set_si64(ibz_t *x, int64_t v)
{
    if (v >= 0) {
        ibz_set_u64(x, (uint64_t)v);
    } else {
        ibz_set_u64(x, (uint64_t)(-v));
        ibz_neg(x, x);
    }
}

void dpe_get_z(ibz_t *x, const dpe_t y)
{
    DPE_EXP_T ey = DPE_EXP(y);

    if (DPE_MANT(y) == 0.0) {
        ibz_init(x);
        return;
    }

    // ey >= DPE_BITSIZE => y는 정수(ULP>=1)
    if (ey >= DPE_BITSIZE) {
        // mant * 2^DPE_BITSIZE 는 (이상적으로) 정수
        double d = (double)(DPE_MANT(y) * DPE_2_POW_BITSIZE);
        uint64_t w = (uint64_t) llround(fabs(d));  // 0 <= w < 2^53

        ibz_set_u64(x, w);

        if (ey > DPE_BITSIZE) {
            ibz_mul_2exp(x, x, (size_t)(ey - DPE_BITSIZE));
        }

        if (DPE_MANT(y) < 0.0) {
            ibz_neg(x, x);
        }
        return;
    }

    // |y| < 1/2 -> 0
    if (ey < 0) {
        ibz_init(x);
        return;
    }

    // 그 외: ldexp 후 반올림 (결과는 <= 2^53 범위)
    double d = (double) DPE_LDEXP(DPE_MANT(y), ey);
    long long v = llround(d);
    ibz_set_si64(x, (int64_t)v);
}

void dpe_set_z(dpe_t x, const ibz_t *y)
{
    if (ibz_is_zero(y)) {
        DPE_MANT(x) = 0.0;
        DPE_EXP(x)  = DPE_EXPMIN;
        return;
    }

    int neg = ibz_is_negative(y);

    ibz_t a;
    ibz_init(&a);
    ibz_abs(&a, y);

    int bits = ibz_bitsize(&a);   // >=1
    int k = (bits < DPE_BITSIZE) ? bits : DPE_BITSIZE; // keep k bits
    int shift = bits - k;         // >=0

    ibz_t top;
    ibz_init(&top);
    ibz_div_2exp(&top, &a, (uint32_t)shift);

    uint64_t hi = top[0];         // top < 2^53 보장
    double mant = (double)hi / ldexp(1.0, k);  // in [0.5,1)
    if (neg) mant = -mant;

    DPE_MANT(x) = (DPE_DOUBLE)mant;
    DPE_EXP(x)  = (DPE_EXP_T)bits;

    ibz_finalize(&a);
    ibz_finalize(&top);
}
