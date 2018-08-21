#include <include/gtest/gtest.h>
#include "../TLwe.h"
#include "../arithmetic.h"
#include <NTL/RR.h>

NTL_CLIENT;


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
        }
    }
}

//void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, UINT64 plaintext_precision) {
TEST(TLWE_TEST, slot_encrypt_decrypt) {
    for (int64_t plaintext_expo : {-5, 10, 20}) {
        for (int64_t level_expo : {20, 63, 128}) {
            for (int64_t N : {300, 400, 500}) {
                for (int64_t plaintext_precision : {10, 20, 40}) {
                    int64_t min_nlimbs = limb_precision(level_expo + plaintext_precision);
                    for (int64_t torus_limbs : {min_nlimbs, min_nlimbs + 1}) {
                        BigTorusParams params(torus_limbs, plaintext_expo, level_expo);
                        TLweParams tlwe_params(N, params);
                        shared_ptr<TLweKey> key = tlwe_keygen(tlwe_params);

                        //fixP load and store to RR
                        for (int trial = 0; trial < 2; trial++) {
                            RR::SetPrecision(torus_limbs * BITS_PER_LIMBS);
                            RR test1 = (random_RR() - 0.5) * power2_RR(plaintext_expo + 1);
                            TLwe c(tlwe_params);
                            slot_encrypt(c, test1, *key, plaintext_precision);

                            RR test2 = slot_decrypt(c, *key);

                            RR::SetOutputPrecision(long(640 * log(2) / log(10)));
                            //cerr << test1 << endl;
                            //cerr << test2 << endl;
                            //cerr << log2Diff(test1, test2)  << " vs. " << plaintext_expo-plaintext_precision << endl;
                            EXPECT_LE(log2Diff(test1, test2), plaintext_expo - plaintext_precision + 1);
                        }
                    }
                }
            }
        }
    }
}
