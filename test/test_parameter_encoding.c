// SPDX-License-Identifier: Apache-2.0

#include <api.h>
#include <encoded_sizes.h>
#include <sig.h>
#include <stdio.h>
#include <string.h>
#include <verification.h>

_Static_assert(CRYPTO_BYTES == SIGNATURE_BYTES, "signature API and encoder sizes disagree");
_Static_assert(CRYPTO_PUBLICKEYBYTES == PUBLICKEY_BYTES, "public key API and encoder sizes disagree");
_Static_assert(CRYPTO_SECRETKEYBYTES == SECRETKEY_BYTES, "secret key API and encoder sizes disagree");

int
main(void)
{
    const unsigned char msg[] = { 42 };
    unsigned char encoded[SIGNATURE_BYTES], roundtrip[SIGNATURE_BYTES], pk_bytes[PUBLICKEY_BYTES];
    signature_t sig = { 0 }, decoded;
    public_key_t pk = { 0 };
    public_key_init(&pk);
    public_key_to_bytes(pk_bytes, &pk);
    public_key_finalize(&pk);
    fp2_set_zero(&sig.E_aux_A);

    sig.backtracking = RESPONSE_COUNTER_BYTES == 2 ? 256 : 42;
    sig.two_resp_length = SQIsign_response_length;
    sig.hint_aux = 11;
    sig.hint_chall = 17;
    sig.chall_coeff[(CHALLENGE_BITS - 1) / RADIX] = (digit_t)1 << ((CHALLENGE_BITS - 1) % RADIX);
    const unsigned int matrix_bit = SQIsign_response_length + 1;
    sig.mat_Bchall_can_to_B_chall[1][1][matrix_bit / RADIX] = (digit_t)1 << (matrix_bit % RADIX);

    signature_to_bytes(encoded, &sig);
    if (!signature_from_bytes(&decoded, encoded) ||
        decoded.backtracking != sig.backtracking || decoded.two_resp_length != sig.two_resp_length ||
        memcmp(decoded.chall_coeff, sig.chall_coeff, sizeof(sig.chall_coeff)) ||
        memcmp(decoded.mat_Bchall_can_to_B_chall, sig.mat_Bchall_can_to_B_chall,
               sizeof(sig.mat_Bchall_can_to_B_chall)) ||
        decoded.hint_aux != sig.hint_aux || decoded.hint_chall != sig.hint_chall) {
        fputs("signature round trip failed\n", stderr);
        return 1;
    }
    if (RESPONSE_COUNTER_BYTES == 2 &&
        (encoded[FP2_ENCODED_BYTES] != 0 || encoded[FP2_ENCODED_BYTES + 1] != 1)) {
        fputs("response counter is not little-endian\n", stderr);
        return 1;
    }
    signature_to_bytes(roundtrip, &decoded);
    if (memcmp(encoded, roundtrip, sizeof(encoded))) {
        fputs("re-encoding changed the signature bytes\n", stderr);
        return 1;
    }

    memset(sig.mat_Bchall_can_to_B_chall, 0, sizeof(sig.mat_Bchall_can_to_B_chall));
    sig.backtracking = SQIsign_response_length + 1;
    sig.two_resp_length = 0;
    signature_to_bytes(encoded, &sig);
    if (sqisign_verify(msg, sizeof(msg), encoded, sizeof(encoded), pk_bytes) != 1) {
        fputs("overlong backtracking was accepted\n", stderr);
        return 1;
    }

    sig.backtracking = 1;
    sig.two_resp_length = SQIsign_response_length;
    signature_to_bytes(encoded, &sig);
    if (sqisign_verify(msg, sizeof(msg), encoded, sizeof(encoded), pk_bytes) != 1) {
        fputs("overlong combined response was accepted\n", stderr);
        return 1;
    }

    unsigned char output = 0;
    unsigned long long mlen = 1;
    if (sqisign_verify(msg, sizeof(msg), encoded, 1, pk_bytes) != 1 ||
        sqisign_open(&output, &mlen, encoded, 1, pk_bytes) != 1 || mlen != 0) {
        fputs("truncated encoding was accepted\n", stderr);
        return 1;
    }

    encoded[FP_FIELD_BITS / 8] |= 1U << (FP_FIELD_BITS % 8);
    if (signature_from_bytes(&decoded, encoded)) {
        fputs("noncanonical field encoding was accepted\n", stderr);
        return 1;
    }

    return 0;
}
