#include <include/gtest/gtest.h>
#include <NTL/RR.h>
#include "../commons.h"
#include "../TRLwe.h"
#include "../arithmetic.h"

NTL_CLIENT;

double log2Diff(const RR &a, const RR &b);

TEST(TRLWE_TEST, trlwe_zero) {
    for (UINT64 nblimbs : {1, 2, 5}) {
        int64_t pubexpo = 1;
        int64_t levexpo = 2;
        for (UINT64 N : {300, 400, 500}) {
            TRLweParams params(N, {nblimbs, pubexpo, levexpo});

            TRLwe c(params);

            zero(c);

            //check zero
            for (UINT64 i = 0; i < N; ++i) {
                for (UINT64 j = 0; j < nblimbs; ++j) {
                    ASSERT_EQ(c.a[0].getAT(i).limbs_end[-1 - j], 0u);
                    ASSERT_EQ(c.a[1].getAT(i).limbs_end[-1 - j], 0u);
                }
            }
        }
    }
}

TEST(TRLWE_TEST, add_noise) {
    for (UINT64 nblimbs : {1, 2, 5}) {
        int64_t pubexpo = 1;
        int64_t levexpo = 2;
        for (UINT64 N : {256, 512, 1024}) {
            TRLweParams params(N, {nblimbs, pubexpo, levexpo});

            BigTorusPolynomial plaintext(N, params.fixp_params);
            BigTorusPolynomial plaintext2(N, params.fixp_params);

            random(plaintext, nblimbs);
            ::copy(plaintext2, plaintext, nblimbs);

            int64_t noise_bits = UINT64(0.75 * nblimbs * BITS_PER_LIMBS);

            add_noise(plaintext2, noise_bits, nblimbs);

            //check randomness of a
            std::set<UINT64> dist;
            for (UINT64 i = 0; i < N; i++) {
                dist.insert(plaintext.getAT(i).limbs_end[-1]);
            }
            ASSERT_GE(dist.size(), 9 * N / 10);

            for (UINT64 i = 0; i < N; i++) {
                RR::SetPrecision(nblimbs * BITS_PER_LIMBS);
                RR p1 = to_RR(plaintext.getAT(i));
                RR p2 = to_RR(plaintext2.getAT(i));
                //cout << p1 << endl;
                //cout << p2 << endl;
                ASSERT_LE(log2Diff(p1, p2), -noise_bits + 5);
                ASSERT_GE(log2Diff(p1, p2), -noise_bits - 15);
            }
        }
    }
}


TEST(TRLWE_TEST, trlwe_encrypt_decrypt_native) {
    for (UINT64 nblimbs : {1, 2, 5}) {
        int64_t pubexpo = 1;
        int64_t levexpo = 2;
        for (UINT64 N : {32, 64, 128}) {
            //cout << "test: " << nblimbs << " " << N << endl;
            TRLweParams params(N, {nblimbs, pubexpo, levexpo});

            std::shared_ptr<TLweKey> key = tlwe_keygen(params);

            TRLwe c(params);
            BigTorusPolynomial plaintext(N, params.fixp_params);
            BigTorusPolynomial plaintext2(N, params.fixp_params);

            random(plaintext, nblimbs);

            int64_t noise_bits = UINT64(0.75 * nblimbs * BITS_PER_LIMBS);
            native_encrypt(c, plaintext, *key, noise_bits);

            //check randomness of a
            std::set<UINT64> dist;
            for (UINT64 i = 0; i < N; i++) {
                dist.insert(c.a[0].getAT(i).limbs_end[-1]);
            }
            ASSERT_GE(dist.size(), 9 * N / 10);

            //check decryption
            native_phase(plaintext2, c, *key, noise_bits);

            for (UINT64 i = 0; i < N; i++) {
                RR::SetPrecision(nblimbs * BITS_PER_LIMBS);
                RR p1 = to_RR(plaintext.getAT(i));
                RR p2 = to_RR(plaintext2.getAT(i));
                //cout << p1 << endl;
                //cout << p2 << endl;
                ASSERT_LE(log2Diff(p1, p2), -noise_bits + 5);
                ASSERT_GE(log2Diff(p1, p2), -noise_bits - 20);
            }
        }
    }
}
