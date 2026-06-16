// SPDX-License-Identifier: Apache-2.0

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <api.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>
#if defined(TARGET_BIG_ENDIAN)
#include <tutil.h>
#endif

typedef enum { PHASE_ALL = 0, PHASE_KEYPAIR, PHASE_SIGN, PHASE_VERIFY } phase_t;

void
bench(size_t runs, phase_t only)
{
    const size_t m_len = 32;
    const size_t sm_len = CRYPTO_BYTES + m_len;

    unsigned char *pkbuf = calloc(runs, CRYPTO_PUBLICKEYBYTES);
    unsigned char *skbuf = calloc(runs, CRYPTO_SECRETKEYBYTES);
    unsigned char *smbuf = calloc(runs, sm_len);
    unsigned char *mbuf = calloc(runs, m_len);

    unsigned char *pk[runs], *sk[runs], *sm[runs], *m[runs];
    for (size_t i = 0; i < runs; ++i) {
        pk[i] = pkbuf + i * CRYPTO_PUBLICKEYBYTES;
        sk[i] = skbuf + i * CRYPTO_SECRETKEYBYTES;
        sm[i] = smbuf + i * sm_len;
        m[i] = mbuf + i * m_len;
        if (randombytes(m[i], m_len))
            abort();
    }

    unsigned long long len;

    printf("%s (%zu iterations, phase=%s)\n", CRYPTO_ALGNAME, runs,
           only == PHASE_KEYPAIR ? "keypair" :
           only == PHASE_SIGN    ? "sign"    :
           only == PHASE_VERIFY  ? "verify"  : "all");

    // Keypair: time it when requested (only=all or only=keypair); otherwise run
    // it untimed because sign/verify need fresh keypairs as input.
    if (only == PHASE_ALL || only == PHASE_KEYPAIR) {
        BENCH_CODE_1(runs);
        crypto_sign_keypair(pk[i], sk[i]);
        BENCH_CODE_2("keypair");
    } else {
        for (size_t i = 0; i < runs; ++i)
            crypto_sign_keypair(pk[i], sk[i]);
    }

    if (only == PHASE_KEYPAIR)
        goto cleanup;

    if (only == PHASE_ALL || only == PHASE_SIGN) {
        BENCH_CODE_1(runs);
        len = sm_len;
        crypto_sign(sm[i], &len, m[i], m_len, sk[i]);
        if (len != sm_len)
            abort();
        BENCH_CODE_2("sign");
    } else {
        for (size_t i = 0; i < runs; ++i) {
            len = sm_len;
            crypto_sign(sm[i], &len, m[i], m_len, sk[i]);
            if (len != sm_len)
                abort();
        }
    }

    if (only == PHASE_SIGN)
        goto cleanup;

    int ret;
    BENCH_CODE_1(runs);
    len = m_len;
    ret = crypto_sign_open(m[i], &len, sm[i], sm_len, pk[i]);
    if (ret)
        abort();
    BENCH_CODE_2("verify");

cleanup:
    free(pkbuf);
    free(skbuf);
    free(smbuf);
    free(mbuf);
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    phase_t only = PHASE_ALL;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }

        if (strcmp(argv[i], "--only=keypair") == 0) { only = PHASE_KEYPAIR; continue; }
        if (strcmp(argv[i], "--only=sign")    == 0) { only = PHASE_SIGN;    continue; }
        if (strcmp(argv[i], "--only=verify")  == 0) { only = PHASE_VERIFY;  continue; }
        if (strcmp(argv[i], "--only=all")     == 0) { only = PHASE_ALL;     continue; }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>] [--only=<phase>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        printf("Where <phase> is one of: keypair, sign, verify, all (default: all).\n"
               "  --only=keypair  : time only keypair (skip sign/verify)\n"
               "  --only=sign     : keypair untimed, time sign, skip verify\n"
               "  --only=verify   : keypair+sign untimed, time only verify\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    cpucycles_init();

    bench(iterations, only);

    return 0;
}
