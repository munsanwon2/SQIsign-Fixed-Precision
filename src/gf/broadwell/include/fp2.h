/* SPDX-License-Identifier: Apache-2.0 */

#ifndef FP2_H
#define FP2_H

#define NO_FP2X_MUL
#define NO_FP2X_SQR

#include <fp2x.h>

void fp2_sq_c0(fp2_t *out, const fp2_t *in);
void fp2_sq_c1(fp_t *out, const fp2_t *in);
void fp2_mul_c0(fp_t *out, const fp2_t *in0, const fp2_t *in1);
void fp2_mul_c1(fp_t *out, const fp2_t *in0, const fp2_t *in1);

static inline void
fp2_mul(fp2_t *x, const fp2_t *y, const fp2_t *z)
{
    fp_t re;
    fp2_mul_c0(&re, y, z);
    fp2_mul_c1(&x->im, y, z);
    x->re = re;
}

static inline void
fp2_sqr(fp2_t *x, const fp2_t *y)
{
    fp2_t tmp;
    fp2_sq_c0(&tmp, y);
    fp2_sq_c1(&x->im, y);
    x->re = tmp.re;
}

#endif
