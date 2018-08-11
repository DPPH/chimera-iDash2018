#include <include/gtest/gtest.h>
#include "../TLwe.h"
#include "../arithmetic.h"
#include <NTL/RR.h>

NTL_CLIENT;

double log2Diff(const RR &a, const RR &b);


TEST(TLWE_TEST, tlwe_keygen) {
    UINT64 nblimbs = 5;
    int64_t pubexpo = 1;
    int64_t levexpo = 2;
    for (UINT64 N : {300, 400, 500}) {
        TLweParams params(N, {nblimbs, pubexpo, levexpo});

        std::shared_ptr<TLweKey> key = tlwe_keygen(params);

        UINT64 count = 0;
        for (UINT64 i = 0; i < N; i++) {
            if (key->key[i]) count++;
        }
        ASSERT_LE(count, 3 * N / 4);
        ASSERT_GE(count, 1 * N / 4);
    }
}

TEST(TLWE_TEST, tlwe_zero) {
    for (UINT64 nblimbs : {1, 2, 5}) {
        int64_t pubexpo = 1;
        int64_t levexpo = 2;
        for (UINT64 N : {300, 400, 500}) {
            TLweParams params(N, {nblimbs, pubexpo, levexpo});

            TLwe c(params);

            zero(c);

            //check zero
            for (UINT64 i = 0; i < N; ++i) {
                for (UINT64 j = 0; j < nblimbs; ++j) {
                    ASSERT_EQ(c.getAT(i).limbs_end[-1 - j], 0u);
                }
            }
        }
    }
}

TEST(TLWE_TEST, tlwe_encrypt_decrypt_native) {
    for (UINT64 nblimbs : {1, 2, 5}) {
        int64_t pubexpo = 1;
        int64_t levexpo = 2;
        for (UINT64 N : {300, 400, 500}) {
            TLweParams params(N, {nblimbs, pubexpo, levexpo});

            std::shared_ptr<TLweKey> key = tlwe_keygen(params);

            TLwe c(params);
            BigTorus plaintext(params.fixp_params);
            BigTorus plaintext2(params.fixp_params);

            random(plaintext);

            int64_t noise_bits = UINT64(0.75 * nblimbs * BITS_PER_LIMBS);
            native_encrypt(c, plaintext, *key, noise_bits);

            //check randomness
            std::set<UINT64> dist;
            for (UINT64 i = 0; i < N; i++) {
                dist.insert(c.getAT(i).limbs_end[-1]);
            }
            ASSERT_GE(dist.size(), 9 * N / 10);

            //check decryption
            native_phase(plaintext2, c, *key);

            RR::SetPrecision(nblimbs * BITS_PER_LIMBS);
            RR p1 = fixp_to_RR(plaintext);
            RR p2 = fixp_to_RR(plaintext2);
            //cout << p1 << endl;
            //cout << p2 << endl;
            ASSERT_LE(log2Diff(p1, p2), -noise_bits + 5);
            ASSERT_GE(log2Diff(p1, p2), -noise_bits - 5);
        }
    }
}

