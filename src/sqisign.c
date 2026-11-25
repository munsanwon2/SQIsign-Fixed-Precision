#include <sig.h>
#include <string.h>
#include <encoded_sizes.h>
#include <verification.h>
#if defined(ENABLE_SIGN)
#include <signature.h>
#endif

#if defined(ENABLE_SIGN)
SQISIGN_API
int
sqisign_keypair(unsigned char *pk, unsigned char *sk)
{
    int ret = 0;
    secret_key_t skt;
    public_key_t pkt = { 0 };
    secret_key_init(&skt);

    ret = !protocols_keygen(&pkt, &skt);

    secret_key_to_bytes(sk, &skt, &pkt);
    public_key_to_bytes(pk, &pkt);
    secret_key_finalize(&skt);
    return ret;
}

SQISIGN_API
int
sqisign_sign(unsigned char *sm,
             unsigned long long *smlen,
             const unsigned char *m,
             unsigned long long mlen,
             const unsigned char *sk)
{
    int ret = 0;
    secret_key_t skt;
    public_key_t pkt = { 0 };
    signature_t sigt;
    secret_key_init(&skt);
    if (!secret_key_from_bytes(&skt, &pkt, sk)) {
        ret = 1;
        *smlen = 0;
        goto err;
    }

    memmove(sm + SIGNATURE_BYTES, m, mlen);

    ret = !protocols_sign(&sigt, &pkt, &skt, sm + SIGNATURE_BYTES, mlen);
    if (ret != 0) {
        *smlen = 0;
        goto err;
    }

    signature_to_bytes(sm, &sigt);
    *smlen = SIGNATURE_BYTES + mlen;

err:
    secret_key_finalize(&skt);
    return ret;
}
#endif

SQISIGN_API
int
sqisign_open(unsigned char *m,
             unsigned long long *mlen,
             const unsigned char *sm,
             unsigned long long smlen,
             const unsigned char *pk)
{
    if (smlen < SIGNATURE_BYTES) {
        *mlen = 0;
        return 1;
    }
    int ret = sqisign_verify(sm + SIGNATURE_BYTES, smlen - SIGNATURE_BYTES, sm, SIGNATURE_BYTES, pk);

    if (!ret) {
        *mlen = smlen - SIGNATURE_BYTES;
        memmove(m, sm + SIGNATURE_BYTES, *mlen);
    } else {
        *mlen = 0;
        memset(m, 0, smlen - SIGNATURE_BYTES);
    }

    return ret;
}

SQISIGN_API
int
sqisign_verify(const unsigned char *m,
               unsigned long long mlen,
               const unsigned char *sig,
               unsigned long long siglen,
               const unsigned char *pk)
{

    if (siglen != SIGNATURE_BYTES)
        return 1;

    public_key_t pkt = { 0 };
    signature_t sigt;

    if (!public_key_from_bytes(&pkt, pk) || !signature_from_bytes(&sigt, sig))
        return 1;

    return !protocols_verify(&sigt, &pkt, m, mlen);
}
