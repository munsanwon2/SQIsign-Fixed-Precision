/*
 * Fixed-precision integer arithmetic implementation
 * Using 2's complement representation for signed integers
 * Supports NIST-I (110), NIST-III (168), NIST-V (222) security levels
 */

#include "intbig.h"

const uint64_t ibz_const_zero[IBZ_LIMBS] = { 0 };
const uint64_t ibz_const_one[IBZ_LIMBS] = { 1 };
const uint64_t ibz_const_two[IBZ_LIMBS] = { 2 };
const uint64_t ibz_const_three[IBZ_LIMBS] = { 3 };

#define KARATSUBA_THRESHOLD 32

// Helper macros
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Bit manipulation helpers
static inline int
clz64(uint64_t x)
{
    if (x == 0)
        return 64;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#else
    int n = 0;
    if (x <= 0x00000000FFFFFFFFULL) {
        n += 32;
        x <<= 32;
    }
    if (x <= 0x0000FFFFFFFFFFFFULL) {
        n += 16;
        x <<= 16;
    }
    if (x <= 0x00FFFFFFFFFFFFFFULL) {
        n += 8;
        x <<= 8;
    }
    if (x <= 0x0FFFFFFFFFFFFFFFULL) {
        n += 4;
        x <<= 4;
    }
    if (x <= 0x3FFFFFFFFFFFFFFFULL) {
        n += 2;
        x <<= 2;
    }
    if (x <= 0x7FFFFFFFFFFFFFFFULL) {
        n += 1;
    }
    return n;
#endif
}

static inline int
ctz64(uint64_t x)
{
    if (x == 0)
        return 64;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    if ((x & 0x00000000FFFFFFFFULL) == 0) {
        n += 32;
        x >>= 32;
    }
    if ((x & 0x000000000000FFFFULL) == 0) {
        n += 16;
        x >>= 16;
    }
    if ((x & 0x00000000000000FFULL) == 0) {
        n += 8;
        x >>= 8;
    }
    if ((x & 0x000000000000000FULL) == 0) {
        n += 4;
        x >>= 4;
    }
    if ((x & 0x0000000000000003ULL) == 0) {
        n += 2;
        x >>= 2;
    }
    if ((x & 0x0000000000000001ULL) == 0) {
        n += 1;
    }
    return n;
#endif
}

// 64x64 -> 128 bit multiplication without __int128
static void
mul64_128(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
    uint32_t a_lo = (uint32_t)a;
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)b;
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64_t p00 = (uint64_t)a_lo * b_lo;
    uint64_t p01 = (uint64_t)a_lo * b_hi;
    uint64_t p10 = (uint64_t)a_hi * b_lo;
    uint64_t p11 = (uint64_t)a_hi * b_hi;

    uint64_t middle = p01 + p10;
    uint64_t carry = (middle < p01) ? (1ULL << 32) : 0;

    *lo = p00 + (middle << 32);
    carry += (*lo < p00) ? 1 : 0;
    *hi = p11 + (middle >> 32) + carry;
}

// Initialize/finalize
void
ibz_init(ibz_t *x)
{
    memset(*x, 0, sizeof(ibz_t));
}

void
ibz_finalize(ibz_t *x)
{
    memset(*x, 0, sizeof(ibz_t));
}

// Copy and swap
void
ibz_copy(ibz_t *target, const ibz_t *value)
{
    memcpy(*target, *value, sizeof(ibz_t));
}

void
ibz_swap(ibz_t *a, ibz_t *b)
{
    uint64_t tmp[IBZ_LIMBS];
    memcpy(tmp, *a, sizeof(ibz_t));
    memcpy(*a, *b, sizeof(ibz_t));
    memcpy(*b, tmp, sizeof(ibz_t));
}

// Check if negative (2's complement)
int
ibz_is_negative(const ibz_t *x)
{
    return ((*x)[IBZ_LIMBS - 1] >> 63) & 1;
}

// Negation (2's complement)
void
ibz_neg(ibz_t *neg, const ibz_t *a)
{
    uint64_t carry = 1;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = ~(*a)[i];
        (*neg)[i] = tmp + carry;
        carry = ((*neg)[i] < tmp) ? 1 : 0;
    }
}

// Absolute value
void
ibz_abs(ibz_t *abs, const ibz_t *a)
{
    if (ibz_is_negative(a)) {
        ibz_neg(abs, a);
    } else {
        ibz_copy(abs, a);
    }
}

// Addition
void
ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b)
{
    uint64_t carry = 0;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = (*a)[i] + carry;
        carry = (tmp < (*a)[i]) ? 1 : 0;
        (*sum)[i] = tmp + (*b)[i];
        carry += ((*sum)[i] < tmp) ? 1 : 0;
    }
}

// Subtraction
void
ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b)
{
    uint64_t borrow = 0;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = (*a)[i] - borrow;
        borrow = (tmp > (*a)[i]) ? 1 : 0;
        uint64_t tmp2 = tmp - (*b)[i];
        borrow += (tmp2 > tmp) ? 1 : 0;
        (*diff)[i] = tmp2;
    }
}

// static void
// ibz_mul_karatsuba_internal(uint64_t *result,
//                            const uint64_t *a,
//                            int a_size,
//                            const uint64_t *b,
//                            int b_size,
//                            uint64_t *workspace)
// {
//     // Base case: use schoolbook for small sizes
//     if (a_size < KARATSUBA_THRESHOLD || b_size < KARATSUBA_THRESHOLD) {
//         // Schoolbook multiplication
//         memset(result, 0, (a_size + b_size) * sizeof(uint64_t));

//         for (int i = 0; i < a_size; i++) {
//             if (a[i] == 0)
//                 continue;

//             uint64_t carry = 0;
//             for (int j = 0; j < b_size; j++) {
//                 uint64_t hi, lo;
//                 mul64_128(a[i], b[j], &hi, &lo);

//                 // Add lo to result[i+j]
//                 result[i + j] += lo;
//                 uint64_t c1 = (result[i + j] < lo) ? 1 : 0;

//                 // Add previous carry
//                 result[i + j] += carry;
//                 c1 += (result[i + j] < carry) ? 1 : 0;

//                 // New carry = hi + c1
//                 carry = hi + c1;
//             }
//             if (carry && (i + b_size < a_size + b_size)) {
//                 result[i + b_size] = carry;
//             }
//         }
//         return;
//     }

//     // Karatsuba: split numbers in half
//     int m = MAX(a_size, b_size) / 2;

//     // a = a1*B^m + a0
//     // b = b1*B^m + b0
//     // result = a1*b1*B^(2m) + ((a1+a0)*(b1+b0) - a1*b1 - a0*b0)*B^m + a0*b0

//     int a0_size = MIN(m, a_size);
//     int a1_size = a_size > m ? a_size - m : 0;
//     int b0_size = MIN(m, b_size);
//     int b1_size = b_size > m ? b_size - m : 0;

//     // Recursive calls for a0*b0 and a1*b1
//     ibz_mul_karatsuba_internal(workspace, a, a0_size, b, b0_size, workspace + 2 * (a_size + b_size));

//     if (a1_size > 0 && b1_size > 0) {
//         ibz_mul_karatsuba_internal(
//             workspace + 2 * m, a + m, a1_size, b + m, b1_size, workspace + 2 * (a_size + b_size));
//     }

//     // Copy results to final location
//     memcpy(result, workspace, (a_size + b_size) * sizeof(uint64_t));
// }

// // OPTIMIZED: Multiplication with Karatsuba for large numbers
// void
// ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b)
// {
//     // Find actual sizes
//     int a_size = IBZ_LIMBS;
//     int b_size = IBZ_LIMBS;

//     while (a_size > 1 && (*a)[a_size - 1] == 0)
//         a_size--;
//     while (b_size > 1 && (*b)[b_size - 1] == 0)
//         b_size--;

//     // // OPTIMIZATION: Use Karatsuba for large multiplications
//     // if (a_size > KARATSUBA_THRESHOLD && b_size > KARATSUBA_THRESHOLD) {
//     //     static uint64_t workspace[6 * IBZ_LIMBS]; // Working space for Karatsuba
//     //     static uint64_t temp_result[2 * IBZ_LIMBS];

//     //     memset(temp_result, 0, sizeof(temp_result));
//     //     ibz_mul_karatsuba_internal(temp_result, *a, a_size, *b, b_size, workspace);

//     //     // Copy result
//     //     for (int i = 0; i < IBZ_LIMBS; i++) {
//     //         (*prod)[i] = temp_result[i];
//     //     }
//     //     return;
//     // }

//     // Fallback to schoolbook for smaller numbers
//     static uint64_t temp_result[2 * IBZ_LIMBS];
//     memset(temp_result, 0, sizeof(temp_result));

//     for (int i = 0; i < a_size; i++) {
//         uint64_t a_limb = (*a)[i];
//         if (a_limb == 0)
//             continue;

//         for (int j = 0; j < b_size; j++) {
//             uint64_t b_limb = (*b)[j];
//             if (b_limb == 0)
//                 continue;

//             if (i + j >= 2 * IBZ_LIMBS)
//                 continue;

//             uint64_t hi, lo;
//             mul64_128(a_limb, b_limb, &hi, &lo);

//             int k = i + j;
//             temp_result[k] += lo;
//             uint64_t carry = (temp_result[k] < lo) ? 1 : 0;

//             k++;
//             if (k < 2 * IBZ_LIMBS) {
//                 temp_result[k] += hi;
//                 uint64_t c2 = (temp_result[k] < hi) ? 1 : 0;
//                 temp_result[k] += carry;
//                 c2 += (temp_result[k] < carry) ? 1 : 0;
//                 carry = c2;

//                 k++;
//                 while (carry && k < 2 * IBZ_LIMBS) {
//                     temp_result[k] += carry;
//                     carry = (temp_result[k] == 0) ? 1 : 0;
//                     k++;
//                 }
//             }
//         }
//     }

//     for (int i = 0; i < IBZ_LIMBS; i++) {
//         (*prod)[i] = temp_result[i];
//     }
// }

void
ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b)
{
    // alias-safe (prod==a or prod==b 대응)
    ibz_t a_in, b_in;
    const ibz_t *A = a, *B = b;
    if (prod == a) { ibz_copy(&a_in, a); A = &a_in; }
    if (prod == b) { ibz_copy(&b_in, b); B = &b_in; }

    int a_neg = ibz_is_negative(A);
    int b_neg = ibz_is_negative(B);
    int neg = (a_neg != b_neg);

    ibz_t aa, bb;
    ibz_abs(&aa, A);
    ibz_abs(&bb, B);

    // Find actual sizes on |a|, |b|
    int a_size = IBZ_LIMBS;
    int b_size = IBZ_LIMBS;

    while (a_size > 1 && aa[a_size - 1] == 0) a_size--;
    while (b_size > 1 && bb[b_size - 1] == 0) b_size--;

    static uint64_t temp_result[2 * IBZ_LIMBS];
    memset(temp_result, 0, sizeof(temp_result));

    for (int i = 0; i < a_size; i++) {
        uint64_t a_limb = aa[i];
        if (a_limb == 0) continue;

        for (int j = 0; j < b_size; j++) {
            uint64_t b_limb = bb[j];
            if (b_limb == 0) continue;

            if (i + j >= 2 * IBZ_LIMBS) continue;

            uint64_t hi, lo;
            mul64_128(a_limb, b_limb, &hi, &lo);

            int k = i + j;
            temp_result[k] += lo;
            uint64_t carry = (temp_result[k] < lo) ? 1 : 0;

            k++;
            if (k < 2 * IBZ_LIMBS) {
                temp_result[k] += hi;
                uint64_t c2 = (temp_result[k] < hi) ? 1 : 0;
                temp_result[k] += carry;
                c2 += (temp_result[k] < carry) ? 1 : 0;
                carry = c2;

                k++;
                while (carry && k < 2 * IBZ_LIMBS) {
                    temp_result[k] += carry;
                    carry = (temp_result[k] == 0) ? 1 : 0;
                    k++;
                }
            }
        }
    }

    // write low limbs
    for (int i = 0; i < IBZ_LIMBS; i++) {
        (*prod)[i] = temp_result[i];
    }

    // apply sign
    if (neg && !ibz_is_zero(prod)) {
        ibz_neg(prod, prod);
    }
}

// void
// ibz_mul_2exp(ibz_t* result, const ibz_t* a, size_t shift)
// {
//     if (shift == 0) {
//         ibz_copy(result, a);
//         return;
//     }

//     size_t limb_shift = shift / 64;
//     size_t bit_shift = shift % 64;

//     ibz_init(result);

//     if (bit_shift == 0) {
//         // limb ���� ����Ʈ��
//         for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
//             (*result)[i + limb_shift] = (*a)[i];
//         }
//     } else {
//         // limb + bit ����Ʈ
//         for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
//             (*result)[i + limb_shift] |= (*a)[i] << bit_shift;
//             if (i + limb_shift + 1 < IBZ_LIMBS) {
//                 (*result)[i + limb_shift + 1] |= (*a)[i] >> (64 - bit_shift);
//             }
//         }
//     }
// }
void
ibz_mul_2exp(ibz_t *result, const ibz_t *a, size_t shift)
{
    // in-place 호출 대비: result == a 이면 따로 복사해서 사용
    ibz_t a_copy;
    const ibz_t *src = a;

    if (result == a) {
        ibz_copy(&a_copy, a);
        src = &a_copy;
    }

    if (shift == 0) {
        ibz_copy(result, src);
        return;
    }

    size_t limb_shift = shift / 64;
    size_t bit_shift  = shift % 64;

    ibz_init(result);

    if (limb_shift >= IBZ_LIMBS) {
        // 전체 범위를 벗어나면 그냥 0
        return;
    }

    if (bit_shift == 0) {
        // limb 단위 shift
        for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
            (*result)[i + limb_shift] = (*src)[i];
        }
    } else {
        // limb + bit shift
        for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
            (*result)[i + limb_shift] |= (*src)[i] << bit_shift;
            if (i + limb_shift + 1 < IBZ_LIMBS) {
                (*result)[i + limb_shift + 1] |= (*src)[i] >> (64 - bit_shift);
            }
        }
    }
}

// Division by power of 2
void
ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp)
{
    if (exp >= IBZ_BITS) {
        if (ibz_is_negative(a)) {
            memset(*quotient, 0xFF, sizeof(ibz_t));
        } else {
            ibz_init(quotient);
        }
        return;
    }

    int limb_shift = exp / 64;
    int bit_shift = exp % 64;

    if (bit_shift == 0) {
        for (int i = 0; i < IBZ_LIMBS - limb_shift; i++) {
            (*quotient)[i] = (*a)[i + limb_shift];
        }
        uint64_t sign = ibz_is_negative(a) ? UINT64_MAX : 0;
        for (int i = IBZ_LIMBS - limb_shift; i < IBZ_LIMBS; i++) {
            (*quotient)[i] = sign;
        }
    } else {
        for (int i = 0; i < IBZ_LIMBS - limb_shift - 1; i++) {
            (*quotient)[i] = ((*a)[i + limb_shift] >> bit_shift) | ((*a)[i + limb_shift + 1] << (64 - bit_shift));
        }
        uint64_t sign = ibz_is_negative(a) ? UINT64_MAX : 0;
        (*quotient)[IBZ_LIMBS - limb_shift - 1] = ((*a)[IBZ_LIMBS - 1] >> bit_shift) | (sign << (64 - bit_shift));
        for (int i = IBZ_LIMBS - limb_shift; i < IBZ_LIMBS; i++) {
            (*quotient)[i] = sign;
        }
    }
}

// Comparison
int
ibz_cmp(const ibz_t *a, const ibz_t *b)
{
    int a_neg = ibz_is_negative(a);
    int b_neg = ibz_is_negative(b);

    if (a_neg != b_neg) {
        return a_neg ? -1 : 1;
    }

    for (int i = IBZ_LIMBS - 1; i >= 0; i--) {
        if ((*a)[i] > (*b)[i])
            return 1;
        if ((*a)[i] < (*b)[i])
            return -1;
    }
    return 0;
}

// Basic predicates
int
ibz_is_zero(const ibz_t *x)
{
    for (int i = 0; i < IBZ_LIMBS; i++) {
        if ((*x)[i] != 0)
            return 0;
    }
    return 1;
}

int
ibz_is_one(const ibz_t *x)
{
    if ((*x)[0] != 1)
        return 0;
    for (int i = 1; i < IBZ_LIMBS; i++) {
        if ((*x)[i] != 0)
            return 0;
    }
    return 1;
}

int
ibz_is_even(const ibz_t *x)
{
    return ((*x)[0] & 1) == 0;
}

int
ibz_is_odd(const ibz_t *x)
{
    return ((*x)[0] & 1) == 1;
}

// Set from int32
void
ibz_set(ibz_t *i, int32_t x)
{
    if (x >= 0) {
        (*i)[0] = (uint64_t)x;
        memset(&(*i)[1], 0, (IBZ_LIMBS - 1) * sizeof(uint64_t));
    } else {
        (*i)[0] = (uint64_t)x;
        memset(&(*i)[1], 0xFF, (IBZ_LIMBS - 1) * sizeof(uint64_t));
    }
}

void
ibz_set_u64(ibz_t *i, uint64_t x)
{
    (*i)[0] = x;
    memset(&(*i)[1], 0, (IBZ_LIMBS - 1) * sizeof(uint64_t));
}

// int
// ibz_convert_to_str(const ibz_t *i, char *str, int base)
// {
//     if (!str || (base != 10 && base != 16))
//         return 0;

//     ibz_t abs_i, base_ibz, q, r;
//     ibz_abs(&abs_i, i);
//     ibz_set(&base_ibz, base);

//     char temp[4096];
//     int pos = 0;

//     if (ibz_is_zero(i)) {
//         str[0] = '0';
//         str[1] = '\0';
//         return 1;
//     }

//     ibz_copy(&q, &abs_i);
//     while (!ibz_is_zero(&q)) {
//         ibz_div(&q, &r, &q, &base_ibz);
//         int digit = (int)r[0];
//         temp[pos++] = (digit < 10) ? ('0' + digit) : ('a' + digit - 10);
//     }

//     int offset = 0;
//     if (ibz_is_negative(i)) {
//         str[offset++] = '-';
//     }

//     for (int ind = 0; ind < pos; ind++) {
//         str[offset + ind] = temp[pos - 1 - ind];
//     }
//     str[offset + pos] = '\0';

//     return 1;
// }
int
ibz_convert_to_str(const ibz_t *i, char *str, int base)
{
    if (!str || (base != 10 && base != 16))
        return 0;

    ibz_t abs_i, base_ibz, q, r;
    ibz_abs(&abs_i, i);
    ibz_set(&base_ibz, base);

    char temp[4096];
    int pos = 0;

    if (ibz_is_zero(i)) {
        str[0] = '0';
        str[1] = '\0';
        return 1;
    }

    ibz_copy(&q, &abs_i);
    while (!ibz_is_zero(&q)) {
        ibz_t q_next;
        ibz_div(&q_next, &r, &q, &base_ibz);   // q_next = q / base, r = q % base

        // 나머지가 혹시라도 base 이상이거나 이상한 값이면 방어
        uint64_t limb0 = r[0];
        int digit = (int)(limb0 % (uint64_t)base);

        temp[pos++] = (digit < 10)
                          ? ('0' + digit)
                          : ('a' + digit - 10);

        ibz_copy(&q, &q_next);
    }

    while (pos > 1 && temp[pos - 1] == '0') {
        pos--;
    }

    int offset = 0;
    if (ibz_is_negative(i)) {
        str[offset++] = '-';
    }

    for (int ind = 0; ind < pos; ind++) {
        str[offset + ind] = temp[pos - 1 - ind];
    }
    str[offset + pos] = '\0';

    return 1;
}


void
ibz_print(const ibz_t *num, int base)
{
    // assert(base == 10 || base == 16);

    int num_size = ibz_size_in_base(num, base);
    char num_str[num_size + 3];
    ibz_convert_to_str(num, num_str, base);
    printf("%s", num_str);
}

int
ibz_set_from_str(ibz_t *i, const char *str, int base)
{
    if (!str || (base != 10 && base != 16))
        return 0;

    ibz_init(i);

    int is_negative = 0;
    int pos = 0;

    if (str[0] == '-') {
        is_negative = 1;
        pos = 1;
    } else if (str[0] == '+') {
        pos = 1;
    }

    ibz_t base_ibz, digit_ibz;
    ibz_set(&base_ibz, base);

    while (str[pos] != '\0') {
        char c = str[pos];
        int digit;

        if (c >= '0' && c <= '9') {
            digit = c - '0';
        } else if (base == 16 && c >= 'a' && c <= 'f') {
            digit = c - 'a' + 10;
        } else if (base == 16 && c >= 'A' && c <= 'F') {
            digit = c - 'A' + 10;
        } else {
            return 0; // Invalid character
        }

        if (digit >= base) {
            return 0; // Invalid digit for base
        }

        // i = i * base + digit
        ibz_mul(i, i, &base_ibz);
        ibz_set(&digit_ibz, digit);
        ibz_add(i, i, &digit_ibz);

        pos++;
    }

    if (is_negative) {
        ibz_neg(i, i);
    }

    return 1;
}

// Get as int32
int32_t
ibz_get(const ibz_t *i)
{
    return (int32_t)((*i)[0] & 0xFFFFFFFF);
}

int
ibz_rand_interval(ibz_t *rand, const ibz_t *a, const ibz_t *b)
{
    ibz_t range;
    ibz_sub(&range, b, a);

    if (ibz_cmp_int32(&range, 0) <= 0) {
        ibz_copy(rand, a);
        return 1;
    }

    int len_bits = ibz_bitsize(&range);
    int len_bytes = (len_bits + 7) / 8;
    int len_limbs = (len_bytes + 7) / 8;

    uint64_t mask = (len_bits % 64 == 0) ? UINT64_MAX : ((1ULL << (len_bits % 64)) - 1);

    // Rejection sampling
    for (int tries = 0; tries < 1000; tries++) {
        ibz_init(rand);

        if (randombytes((unsigned char *)(*rand), len_bytes) != 0) {
            return 0;
        }

        if (len_limbs > 0 && len_limbs <= IBZ_LIMBS) {
            (*rand)[len_limbs - 1] &= mask;
        }

        if (ibz_cmp(rand, &range) < 0) {
            ibz_add(rand, rand, a);
            return 1;
        }
    }

    return 0;
}

int
ibz_rand_interval_i(ibz_t *rand, int32_t a, int32_t b)
{
    ibz_t a_ibz, b_ibz;
    ibz_set(&a_ibz, a);
    ibz_set(&b_ibz, b);
    return ibz_rand_interval(rand, &a_ibz, &b_ibz);
}

int
ibz_rand_interval_minm_m(ibz_t *rand, int32_t m)
{
    ibz_t m_big, neg_m;
    ibz_set(&m_big, m);
    ibz_neg(&neg_m, &m_big);
    return ibz_rand_interval(rand, &neg_m, &m_big);
}

int
ibz_rand_interval_bits(ibz_t *rand, uint32_t m)
{
    ibz_t max_val, min_val;
    ibz_set(&max_val, 1);
    ibz_mul_2exp(&max_val, &max_val, m);
    ibz_neg(&min_val, &max_val);
    return ibz_rand_interval(rand, &min_val, &max_val);
}

// Compare with int32
int
ibz_cmp_int32(const ibz_t *x, int32_t y)
{
    ibz_t tmp;
    ibz_set(&tmp, y);
    return ibz_cmp(x, &tmp);
}

// Get bit size (absolute value)
int
ibz_bitsize(const ibz_t *a)
{
    // Handle negative numbers by checking absolute value
    if (ibz_is_negative(a)) {
        ibz_t abs_a;
        ibz_neg(&abs_a, a);
        for (int i = IBZ_LIMBS - 1; i >= 0; i--) {
            if (abs_a[i] != 0) {
                return i * 64 + (64 - clz64(abs_a[i]));
            }
        }
        return 0;
    }

    // Positive numbers
    for (int i = IBZ_LIMBS - 1; i >= 0; i--) {
        if ((*a)[i] != 0) {
            return i * 64 + (64 - clz64((*a)[i]));
        }
    }
    return 0;
}

// Get size in base
size_t
ibz_size_in_base(const ibz_t *a, int base)
{
    if (ibz_is_zero(a)) {
        return 1;
    }

    ibz_t abs_a, temp, base_ibz;
    ibz_abs(&abs_a, a);
    ibz_copy(&temp, &abs_a);
    ibz_set(&base_ibz, base);

    size_t count = 0;

    while (!ibz_is_zero(&temp)) {
        ibz_t q, r;
        ibz_div(&q, &r, &temp, &base_ibz);
        ibz_copy(&temp, &q);
        count++;
    }

    return count > 0 ? count : 1;
}

// Convert to digit array
void
ibz_to_digits(digit_t *target, const ibz_t *a)
{
    int actual_limbs = IBZ_LIMBS;
    for (int i = IBZ_LIMBS - 1; i >= 0; i--) {
        if ((*a)[i] != 0) {
            actual_limbs = i + 1;
            break;
        }
    }

    if (actual_limbs == 0 || ibz_is_zero(a)) {
        target[0] = 0;
        return;
    }

    for (int i = 0; i < actual_limbs; i++) {
        target[i] = (*a)[i];
    }
}

// Copy from digit array
void
ibz_copy_digits(ibz_t *a, const digit_t *digits, size_t len)
{
    ibz_init(a);
    for (size_t i = 0; i < len && i < IBZ_LIMBS; i++) {
        (*a)[i] = digits[i];
    }
}

// Get 2-adic valuation (trailing zeros)
int
ibz_two_adic(ibz_t *pow)
{
    for (int i = 0; i < IBZ_LIMBS; i++) {
        if ((*pow)[i] != 0) {
            return i * 64 + ctz64((*pow)[i]);
        }
    }
    return IBZ_BITS;
}

// void
// ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b)
// {
//     if (ibz_is_zero(b)) {
//         ibz_init(quotient);
//         ibz_init(remainder);
//         return;
//     }

//     ibz_init(quotient);
//     ibz_init(remainder);

//     int a_neg = ibz_is_negative(a);
//     int b_neg = ibz_is_negative(b);
//     int quot_neg = (a_neg != b_neg);

//     ibz_t dividend, divisor;
//     ibz_abs(&dividend, a);
//     ibz_abs(&divisor, b);

//     if (ibz_cmp(&dividend, &divisor) < 0) {
//         ibz_copy(remainder, a);
//         return;
//     }

//     // OPTIMIZATION: Use word-by-word division for better performance
//     ibz_t q, r;
//     ibz_init(&q);
//     ibz_copy(&r, &dividend);

//     int divisor_bits = ibz_bitsize(&divisor);
//     int dividend_bits = ibz_bitsize(&dividend);

//     // Shift divisor to align with dividend
//     int shift = dividend_bits - divisor_bits;
//     ibz_t shifted_divisor;
//     ibz_mul_2exp(&shifted_divisor, &divisor, shift);

//     // Division loop
//     for (int i = shift; i >= 0; i--) {
//         if (ibz_cmp(&r, &shifted_divisor) >= 0) {
//             ibz_sub(&r, &r, &shifted_divisor);
//             q[i / 64] |= (1ULL << (i % 64));
//         }
//         if (i > 0) {
//             ibz_div_2exp(&shifted_divisor, &shifted_divisor, 1);
//         }
//     }

//     if (quot_neg && !ibz_is_zero(&q)) {
//         ibz_neg(&q, &q);
//     }

//     if (a_neg && !ibz_is_zero(&r)) {
//         ibz_neg(&r, &r);
//     }

//     ibz_copy(quotient, &q);
//     ibz_copy(remainder, &r);
// }

void
ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b)
{
    ibz_t q, r;
    ibz_t dividend, divisor;

    if (ibz_is_zero(b)) {
        if (quotient)  ibz_init(quotient);
        if (remainder) ibz_init(remainder);
        return;
    }

    int a_neg   = ibz_is_negative(a);
    int b_neg   = ibz_is_negative(b);
    int quot_neg = (a_neg != b_neg);

    ibz_abs(&dividend, a);
    ibz_abs(&divisor,  b);

    // |a| < |b| 인 경우: 몫=0, 나머지=a 그대로
    if (ibz_cmp(&dividend, &divisor) < 0) {
        if (quotient)  ibz_init(quotient);
        if (remainder) ibz_copy(remainder, a);
        return;
    }

    ibz_init(&q);
    ibz_copy(&r, &dividend);

    int divisor_bits  = ibz_bitsize(&divisor);
    int dividend_bits = ibz_bitsize(&dividend);

    int shift = dividend_bits - divisor_bits;
    ibz_t shifted_divisor;
    ibz_mul_2exp(&shifted_divisor, &divisor, shift);

    // 메인 division 루프
    for (int i = shift; i >= 0; i--) {
        if (ibz_cmp(&r, &shifted_divisor) >= 0) {
            ibz_sub(&r, &r, &shifted_divisor);
            q[i / 64] |= (1ULL << (i % 64));
        }
        if (i > 0) {
            ibz_div_2exp(&shifted_divisor, &shifted_divisor, 1);
        }
    }

    // 부호 조정
    if (quot_neg && !ibz_is_zero(&q)) {
        ibz_neg(&q, &q);
    }
    if (a_neg && !ibz_is_zero(&r)) {
        ibz_neg(&r, &r);
    }

    // 마지막에만 결과 내보내기 (alias-safe)
    if (quotient)  ibz_copy(quotient,  &q);
    if (remainder) ibz_copy(remainder, &r);
}


// Floor division
void
ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d)
{
    ibz_div(q, r, n, d);

    if (ibz_is_negative(r) && !ibz_is_zero(r)) {
        if (ibz_is_negative(d)) {
            ibz_t one;
            ibz_set(&one, 1);
            ibz_add(q, q, &one);
            ibz_sub(r, r, d);
        } else {
            ibz_t one;
            ibz_set(&one, 1);
            ibz_sub(q, q, &one);
            ibz_add(r, r, d);
        }
    }
}

// Modulo
void
ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b)
{
    ibz_t q;
    ibz_init(&q);
    ibz_div_floor(&q, r, a, b);
}

// Power
void
ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e)
{
    ibz_t result, base;
    ibz_set(&result, 1);
    ibz_copy(&base, x);

    while (e > 0) {
        if (e & 1) {
            ibz_mul(&result, &result, &base);
        }
        ibz_mul(&base, &base, &base);
        e >>= 1;
    }

    ibz_copy(pow, &result);
}

// Modulo for unsigned long
unsigned long
ibz_mod_ui(const ibz_t *n, unsigned long d)
{
    if (d == 0)
        return 0;

    ibz_t divisor, remainder, quotient;
    ibz_set_u64(&divisor, d);
    ibz_div(&quotient, &remainder, n, &divisor);

    if (ibz_is_negative(&remainder)) {
        ibz_t d_ibz;
        ibz_set_u64(&d_ibz, d);
        ibz_add(&remainder, &remainder, &d_ibz);
    }

    return (unsigned long)remainder[0];
}

// Probabilistic primality test (Miller-Rabin)
int
ibz_probab_prime(const ibz_t *n, int reps)
{
    if (ibz_cmp_int32(n, 2) == 0)
        return 1;
    if (ibz_cmp_int32(n, 3) == 0)
        return 1;
    if (ibz_is_even(n))
        return 0;
    if (ibz_cmp_int32(n, 1) <= 0)
        return 0;

    ibz_t n_minus_1, d, a, x, temp;
    ibz_sub(&n_minus_1, n, (const ibz_t *)&ibz_const_one);
    ibz_copy(&d, &n_minus_1);

    while (ibz_is_even(&d)) {
        ibz_div_2exp(&d, &d, 1);
    }

    for (int i = 0; i < reps; i++) {
        ibz_t two, n_minus_2;
        ibz_set(&two, 2);
        ibz_sub(&n_minus_2, n, &two);
        ibz_rand_interval(&a, &two, &n_minus_2);

        // x = a^d mod n
        ibz_pow_mod(&x, &a, &d, n);

        if (ibz_is_one(&x) || ibz_cmp(&x, &n_minus_1) == 0) {
            continue;
        }

        ibz_copy(&temp, &d);
        int composite = 1;
        while (ibz_cmp(&temp, &n_minus_1) < 0) {
            ibz_mul(&x, &x, &x);
            ibz_mod(&x, &x, n);

            if (ibz_is_one(&x)) {
                return 0; // Composite
            }
            if (ibz_cmp(&x, &n_minus_1) == 0) {
                composite = 0;
                break;
            }

            ibz_mul(&temp, &temp, (const ibz_t *)&ibz_const_two);
        }

        if (composite) {
            return 0;
        }
    }

    return 1; // Probably prime
}

// Modular exponentiation
void
ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m)
{
    ibz_t result, base, exp, temp;
    ibz_set(&result, 1);
    ibz_copy(&base, x);
    ibz_copy(&exp, e);

    ibz_mod(&base, &base, m);

    while (!ibz_is_zero(&exp)) {
        if (ibz_is_odd(&exp)) {
            ibz_mul(&temp, &result, &base);
            ibz_mod(&result, &temp, m);
        }
        ibz_mul(&temp, &base, &base);
        ibz_mod(&base, &temp, m);
        ibz_div_2exp(&exp, &exp, 1);
    }

    ibz_copy(pow, &result);
}

// GCD using binary algorithm
void
ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b)
{
    ibz_t u, v;
    ibz_abs(&u, a);
    ibz_abs(&v, b);

    if (ibz_is_zero(&u)) {
        ibz_copy(gcd, &v);
        return;
    }
    if (ibz_is_zero(&v)) {
        ibz_copy(gcd, &u);
        return;
    }

    int shift = MIN(ibz_two_adic(&u), ibz_two_adic(&v));
    ibz_div_2exp(&u, &u, shift);
    ibz_div_2exp(&v, &v, shift);

    while (!ibz_is_zero(&u)) {
        ibz_div_2exp(&u, &u, ibz_two_adic(&u));
        ibz_div_2exp(&v, &v, ibz_two_adic(&v));

        if (ibz_cmp(&u, &v) > 0) {
            ibz_swap(&u, &v);
        }

        ibz_sub(&v, &v, &u);
    }

    // Restore common factors of 2
    ibz_mul_2exp(&v, &v, shift);

    ibz_copy(gcd, &v);
}

void
ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y,
           const ibz_t *a, const ibz_t *b)
{
    // 특수 케이스: a == 0
    if (ibz_is_zero(a)) {
        ibz_abs(gcd, b);
        ibz_set(x, 0);
        if (ibz_is_negative(b)) {
            ibz_set(y, -1);
        } else {
            ibz_set(y, 1);
        }
        return;
    }

    // 특수 케이스: b == 0
    if (ibz_is_zero(b)) {
        ibz_abs(gcd, a);
        if (ibz_is_negative(a)) {
            ibz_set(x, -1);
        } else {
            ibz_set(x, 1);
        }
        ibz_set(y, 0);
        return;
    }

    ibz_t aa, bb;
    ibz_t x0, x1, y0, y1;
    ibz_t q, r;
    ibz_t tmp1, tmp2;

    // a, b 는 나중에 다시 써야 할 수 있으니 로컬에 복사
    ibz_copy(&aa, a);
    ibz_copy(&bb, b);

    // (x0, y0) = (1, 0), (x1, y1) = (0, 1)
    ibz_set(&x0, 1);
    ibz_set(&y0, 0);
    ibz_set(&x1, 0);
    ibz_set(&y1, 1);

    // 표준 확장 유클리드 알고리즘
    while (!ibz_is_zero(&bb)) {
        // aa = q * bb + r  (ibz_div: aa, bb의 부호에 맞는 q, r)
        ibz_div(&q, &r, &aa, &bb);

        // (aa, bb) <- (bb, r)
        ibz_copy(&aa, &bb);
        ibz_copy(&bb, &r);

        // (x0, x1) <- (x1, x0 - q*x1)
        ibz_mul(&tmp1, &q, &x1);    // tmp1 = q * x1
        ibz_sub(&tmp2, &x0, &tmp1); // tmp2 = x0 - q*x1
        ibz_copy(&x0, &x1);
        ibz_copy(&x1, &tmp2);

        // (y0, y1) <- (y1, y0 - q*y1)
        ibz_mul(&tmp1, &q, &y1);    // tmp1 = q * y1
        ibz_sub(&tmp2, &y0, &tmp1); // tmp2 = y0 - q*y1
        ibz_copy(&y0, &y1);
        ibz_copy(&y1, &tmp2);
    }

    // 이제 aa 가 gcd(a,b) (부호 포함)
    ibz_copy(gcd, &aa);
    ibz_copy(x, &x0);
    ibz_copy(y, &y0);

    // gcd 를 양수로 정규화: gcd < 0 이면 gcd,x,y 모두 부호 반전
    if (ibz_is_negative(gcd)) {
        ibz_neg(gcd, gcd);
        ibz_neg(x, x);
        ibz_neg(y, y);
    }
}


// void
// ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y, const ibz_t *a, const ibz_t *b)
// {
//     if (ibz_is_zero(a)) {
//         ibz_abs(gcd, b);
//         ibz_set(x, 0);
//         if (ibz_is_negative(b)) {
//             ibz_set(y, -1);
//         } else {
//             ibz_set(y, 1);
//         }
//         return;
//     }

//     if (ibz_is_zero(b)) {
//         ibz_abs(gcd, a);
//         if (ibz_is_negative(a)) {
//             ibz_set(x, -1);
//         } else {
//             ibz_set(x, 1);
//         }
//         ibz_set(y, 0);
//         return;
//     }

//     int a_was_neg = ibz_is_negative(a);
//     int b_was_neg = ibz_is_negative(b);

//     ibz_t u, v, A, B, C, D, temp;

//     ibz_abs(&u, a);
//     ibz_abs(&v, b);

//     ibz_set(&A, 1);
//     ibz_set(&B, 0);
//     ibz_set(&C, 0);
//     ibz_set(&D, 1);

//     int shift = 0;
//     while (ibz_is_even(&u) && ibz_is_even(&v)) {
//         ibz_div_2exp(&u, &u, 1);
//         ibz_div_2exp(&v, &v, 1);
//         shift++;
//     }

//     ibz_t orig_u, orig_v;
//     ibz_copy(&orig_u, &u);
//     ibz_copy(&orig_v, &v);

//     // Binary Extended GCD
//     while (!ibz_is_zero(&u)) {
//         while (ibz_is_even(&u)) {
//             ibz_div_2exp(&u, &u, 1);

//             if (ibz_is_even(&A) && ibz_is_even(&B)) {
//                 ibz_div_2exp(&A, &A, 1);
//                 ibz_div_2exp(&B, &B, 1);
//             } else {
//                 // A = (A + orig_v) / 2
//                 ibz_add(&temp, &A, &orig_v);
//                 ibz_div_2exp(&A, &temp, 1);
//                 // B = (B - orig_u) / 2
//                 ibz_sub(&temp, &B, &orig_u);
//                 ibz_div_2exp(&B, &temp, 1);
//             }
//         }

//         while (ibz_is_even(&v)) {
//             ibz_div_2exp(&v, &v, 1);

//             if (ibz_is_even(&C) && ibz_is_even(&D)) {
//                 ibz_div_2exp(&C, &C, 1);
//                 ibz_div_2exp(&D, &D, 1);
//             } else {
//                 // C = (C + orig_v) / 2
//                 ibz_add(&temp, &C, &orig_v);
//                 ibz_div_2exp(&C, &temp, 1);
//                 // D = (D - orig_u) / 2
//                 ibz_sub(&temp, &D, &orig_u);
//                 ibz_div_2exp(&D, &temp, 1);
//             }
//         }

//         if (ibz_cmp(&u, &v) >= 0) {
//             ibz_sub(&u, &u, &v);
//             ibz_sub(&A, &A, &C);
//             ibz_sub(&B, &B, &D);
//         } else {
//             // v = v - u, C = C - A, D = D - B
//             ibz_sub(&v, &v, &u);
//             ibz_sub(&C, &C, &A);
//             ibz_sub(&D, &D, &B);
//         }
//     }

//     // gcd = v * 2^shift
//     ibz_copy(gcd, &v);
//     ibz_mul_2exp(gcd, gcd, shift);

//     ibz_copy(x, &C);
//     ibz_copy(y, &D);

//     if (a_was_neg) {
//         ibz_neg(x, x);
//     }
//     if (b_was_neg) {
//         ibz_neg(y, y);
//     }
// }

// Modular inverse: a^(-1) mod m
int
ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod)
{
    if (ibz_is_zero(mod) || ibz_cmp_int32(mod, 0) <= 0) {
        return 0;
    }

    if (ibz_is_zero(a)) {
        return 0;
    }

    ibz_t gcd, x, y;
    ibz_t a_mod;

    ibz_mod(&a_mod, a, mod);

    if (ibz_is_zero(&a_mod)) {
        return 0;
    }

    ibz_gcdext(&gcd, &x, &y, &a_mod, mod);

    if (!ibz_is_one(&gcd)) {
        return 0;
    }

    ibz_mod(inv, &x, mod);

    return 1;
}

// Check divisibility
int
ibz_divides(const ibz_t *a, const ibz_t *b)
{
    ibz_t q, r;
    ibz_div(&q, &r, a, b);
    return ibz_is_zero(&r);
}

// Integer square root
int
ibz_sqrt(ibz_t *sqrt, const ibz_t *a)
{
    if (ibz_is_negative(a)) {
        return 0;
    }

    if (ibz_is_zero(a)) {
        ibz_init(sqrt);
        return 1;
    }

    ibz_t x, x_prev, temp;

    ibz_copy(&x, a);

    // Newton-Raphson: x_new = (x + a/x) / 2
    for (int iter = 0; iter < 1000; iter++) {
        ibz_copy(&x_prev, &x);

        ibz_t q, r;
        ibz_div(&q, &r, a, &x);
        ibz_add(&temp, &x, &q);
        ibz_div_2exp(&x, &temp, 1);

        if (ibz_cmp(&x, &x_prev) == 0) {
            break;
        }

        if (ibz_cmp(&x, &x_prev) > 0) {
            ibz_copy(&x, &x_prev);
            break;
        }
    }

    ibz_mul(&temp, &x, &x);
    if (ibz_cmp(&temp, a) == 0) {
        ibz_copy(sqrt, &x);
        return 1;
    }

    return 0;
}

// Square root floor
void
ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a)
{
    if (ibz_is_negative(a) || ibz_is_zero(a)) {
        ibz_init(sqrt);
        return;
    }

    ibz_t x, x_next, temp, two;
    ibz_set(&two, 2);

    int bits = ibz_bitsize(a);
    ibz_set(&x, 1);
    ibz_mul_2exp(&x, &x, (bits + 1) / 2);

    // Newton-Raphson with proper termination
    for (int iter = 0; iter < 1000; iter++) {
        ibz_t q, r;
        ibz_div(&q, &r, a, &x);
        ibz_add(&temp, &x, &q);
        ibz_div_2exp(&x_next, &temp, 1);

        if (ibz_cmp(&x_next, &x) >= 0) {
            break;
        }

        ibz_copy(&x, &x_next);
    }

    ibz_copy(sqrt, &x);
}

// Legendre symbol
int
ibz_legendre(const ibz_t *a, const ibz_t *p)
{
    ibz_t a_mod, exp, result, p_minus_1;

    ibz_mod(&a_mod, a, p);

    if (ibz_is_zero(&a_mod)) {
        return 0;
    }

    ibz_sub(&p_minus_1, p, (const ibz_t *)&ibz_const_one);
    ibz_div_2exp(&exp, &p_minus_1, 1);

    ibz_pow_mod(&result, &a_mod, &exp, p);

    if (ibz_is_one(&result)) {
        return 1;
    } else {
        return -1;
    }
}

// Modular square root
// int
// ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p)
// {
//     ibz_t a_mod;
//     ibz_mod(&a_mod, a, p);

//     if (ibz_is_zero(&a_mod)) {
//         ibz_init(sqrt);
//         return 1;
//     }

//     if (ibz_legendre(&a_mod, p) != 1) {
//         return 0;
//     }

//     ibz_t p_mod_4, three, four;
//     ibz_set(&three, 3);
//     ibz_set(&four, 4);
//     ibz_mod(&p_mod_4, p, &four);

//     if (ibz_cmp(&p_mod_4, &three) == 0) {
//         ibz_t exp;
//         ibz_add(&exp, p, (const ibz_t *)&ibz_const_one);
//         ibz_div_2exp(&exp, &exp, 2);
//         ibz_pow_mod(sqrt, &a_mod, &exp, p);
//         return 1;
//     }

//     return 0;
// }

// Modular square root for odd prime p (Tonelli–Shanks)
// returns 1 if sqrt exists (and writes it), 0 otherwise.
int
ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p)
{
    // handle p == 2 (not really used in your setting, but safe)
    if (ibz_cmp_int32(p, 2) == 0) {
        ibz_t am;
        ibz_mod(&am, a, p);
        ibz_copy(sqrt, &am);
        return 1;
    }

    // a := a mod p
    ibz_t n;
    ibz_mod(&n, a, p);
    if (ibz_is_zero(&n)) {
        ibz_init(sqrt);
        return 1;
    }

    // Check quadratic residue via Legendre
    if (ibz_legendre(&n, p) != 1) {
        return 0;
    }

    // Factor p-1 = q * 2^s with q odd
    ibz_t pm1, q;
    ibz_sub(&pm1, p, (const ibz_t *)&ibz_const_one);
    ibz_copy(&q, &pm1);

    int s = 0;
    while (ibz_is_even(&q)) {
        ibz_div_2exp(&q, &q, 1);
        s++;
    }

    // If s == 1, p % 4 == 3 fast-path works too, but TS handles it fine.

    // Find z, a quadratic non-residue mod p
    ibz_t z;
    ibz_set(&z, 2);
    while (ibz_legendre(&z, p) != -1) {
        ibz_add(&z, &z, (const ibz_t *)&ibz_const_one);
    }

    // c = z^q mod p
    ibz_t c;
    ibz_pow_mod(&c, &z, &q, p);

    // x = n^((q+1)/2) mod p
    ibz_t q1, e;
    ibz_add(&q1, &q, (const ibz_t *)&ibz_const_one);
    ibz_div_2exp(&e, &q1, 1);
    ibz_t x;
    ibz_pow_mod(&x, &n, &e, p);

    // t = n^q mod p
    ibz_t t;
    ibz_pow_mod(&t, &n, &q, p);

    int m = s;

    ibz_t one;
    ibz_set(&one, 1);

    while (!ibz_is_one(&t)) {
        // Find least i (0 < i < m) such that t^(2^i) == 1
        int i = 0;
        ibz_t t2i;
        ibz_copy(&t2i, &t);

        for (i = 1; i < m; i++) {
            // t2i = t2i^2 mod p
            ibz_mul(&t2i, &t2i, &t2i);
            ibz_mod(&t2i, &t2i, p);
            if (ibz_is_one(&t2i))
                break;
        }

        // Should always find such i
        if (i == m) {
            return 0; // defensive (shouldn't happen if p prime and residue)
        }

        // b = c^(2^(m-i-1)) mod p
        ibz_t b;
        ibz_copy(&b, &c);
        for (int j = 0; j < (m - i - 1); j++) {
            ibz_mul(&b, &b, &b);
            ibz_mod(&b, &b, p);
        }

        // x = x*b mod p
        ibz_mul(&x, &x, &b);
        ibz_mod(&x, &x, p);

        // t = t*b^2 mod p
        ibz_t b2;
        ibz_mul(&b2, &b, &b);
        ibz_mod(&b2, &b2, p);

        ibz_mul(&t, &t, &b2);
        ibz_mod(&t, &t, p);

        // c = b^2 mod p
        ibz_copy(&c, &b2);

        m = i;
    }

    ibz_copy(sqrt, &x);
    return 1;
}
