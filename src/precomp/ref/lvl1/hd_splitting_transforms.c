#include <hd_splitting_transforms.h>

#define FP2_ZERO 0
#define FP2_ONE 1
#define FP2_I 2
#define FP2_MINUS_ONE 3
#define FP2_MINUS_I 4

const int EVEN_INDEX[10][2] = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 0}, {1, 2}, {2, 0}, {2, 1}, {3, 0}, {3, 3}};
const int CHI_EVAL[4][4] = {{1, 1, 1, 1}, {1, -1, 1, -1}, {1, 1, -1, -1}, {1, -1, -1, 1}};
const fp2_t FP2_CONSTANTS[5] = {{
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
,
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
}, {
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x555555555555555, 0x0, 0x0, 0x0, 0x0, 0x10}}
#elif RADIX == 64
{0x15, 0x0, 0x0, 0x0, 0x0, 0x2000000000000}
#else
{0x555, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x10000}
#endif
,
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
}, {
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
,
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x555555555555555, 0x0, 0x0, 0x0, 0x0, 0x10}}
#elif RADIX == 64
{0x15, 0x0, 0x0, 0x0, 0x0, 0x2000000000000}
#else
{0x555, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x10000}
#endif
}, {
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0xfaaaaaaaaaaaaaaa, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x1f}}
#elif RADIX == 64
{0x7fffffffffffea, 0x7fffffffffffff, 0x7fffffffffffff, 0x7fffffffffffff, 0x7fffffffffffff, 0x3ffffffffffff}
#else
{0xffffaaa, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0x1ffff}
#endif
,
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
}, {
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}
#elif RADIX == 64
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#else
{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}
#endif
,
#if defined(SQISIGN_GF_IMPL_BROADWELL)
{{0xfaaaaaaaaaaaaaaa, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x1f}}
#elif RADIX == 64
{0x7fffffffffffea, 0x7fffffffffffff, 0x7fffffffffffff, 0x7fffffffffffff, 0x7fffffffffffff, 0x3ffffffffffff}
#else
{0xffffaaa, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0x1ffff}
#endif
}};
const precomp_basis_change_matrix_t SPLITTING_TRANSFORMS[10] = {{{{FP2_ONE, FP2_I, FP2_ONE, FP2_I}, {FP2_ONE, FP2_MINUS_I, FP2_MINUS_ONE, FP2_I}, {FP2_ONE, FP2_I, FP2_MINUS_ONE, FP2_MINUS_I}, {FP2_MINUS_ONE, FP2_I, FP2_MINUS_ONE, FP2_I}}}, {{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}, {FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO, FP2_ZERO}}}, {{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}, {FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO}}}, {{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE}}}, {{{FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE}, {FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE}, {FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE}}}, {{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}}}, {{{FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE}, {FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE}}}, {{{FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE}, {FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE}}}, {{{FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE}, {FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE}, {FP2_MINUS_ONE, FP2_ONE, FP2_ONE, FP2_MINUS_ONE}}}, {{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}}}};
const precomp_basis_change_matrix_t NORMALIZATION_TRANSFORMS[6] = {{{{FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}, {FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}}}, {{{FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE}, {FP2_ZERO, FP2_ZERO, FP2_ONE, FP2_ZERO}, {FP2_ZERO, FP2_ONE, FP2_ZERO, FP2_ZERO}, {FP2_ONE, FP2_ZERO, FP2_ZERO, FP2_ZERO}}}, {{{FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE}, {FP2_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE}, {FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE}}}, {{{FP2_ONE, FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE}, {FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE}, {FP2_MINUS_ONE, FP2_ONE, FP2_MINUS_ONE, FP2_ONE}, {FP2_ONE, FP2_ONE, FP2_ONE, FP2_ONE}}}, {{{FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE}, {FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I}, {FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I}, {FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE}}}, {{{FP2_ONE, FP2_I, FP2_I, FP2_MINUS_ONE}, {FP2_I, FP2_ONE, FP2_MINUS_ONE, FP2_I}, {FP2_I, FP2_MINUS_ONE, FP2_ONE, FP2_I}, {FP2_MINUS_ONE, FP2_I, FP2_I, FP2_ONE}}}};
