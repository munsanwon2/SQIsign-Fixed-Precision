#include "hnf_internal.h"
#include "intbig.h"

void
ibz_xgcd(ibz_t *gcd, ibz_t *u, ibz_t *v, const ibz_t *a, const ibz_t *b)
{
    ibz_gcdext(gcd, u, v, a, b);
}