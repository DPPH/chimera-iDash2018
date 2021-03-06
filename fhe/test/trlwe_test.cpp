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

TEST(TRLWE_TEST, fixp_encrypt_decrypt) {
    for (int64_t plaintext_expo : {-5, 10, 20}) {
        for (int64_t level_expo : {20, 63, 128}) {
            for (int64_t N : {32, 64, 256}) {
                for (int64_t plaintext_precision : {10, 20, 40}) {
                    int64_t min_nlimbs = limb_precision(level_expo + plaintext_precision);
                    for (int64_t torus_limbs : {min_nlimbs, min_nlimbs + 1}) {
                        BigTorusParams params(torus_limbs, plaintext_expo, level_expo);
                        TRLweParams trlwe_params(N, params);
                        shared_ptr<TLweKey> key = tlwe_keygen(trlwe_params);

                        //fixP load and store to RR
                        for (int trial = 0; trial < 1; trial++) {
                            RR::SetPrecision(torus_limbs * BITS_PER_LIMBS);
                            vec_RR test1(INIT_SIZE, N);
                            for (int i = 0; i < N; i++) {
                                test1[i] = (random_RR() - 0.5) * power2_RR(plaintext_expo + 1);
                            }
                            TRLwe c(trlwe_params);
                            fixp_encrypt(c, test1, *key, plaintext_precision);

                            vec_RR test2 = fixp_decrypt(c, *key);

                            RR::SetOutputPrecision(long(640 * log(2) / log(10)));
                            //cerr << test1 << endl;
                            //cerr << test2 << endl;
                            for (int i = 0; i < N; i++) {
                                //cerr << log2Diff(test1[i], test2[i])  << " vs. " << plaintext_expo-plaintext_precision << endl;
                                EXPECT_LE(log2Diff(test1[i], test2[i]), plaintext_expo - plaintext_precision + 1);
                            }
                        }
                    }
                }
            }
        }
    }
}


TEST(TRLWE_TEST, pubKSKeyGen) {

}

TEST(TRLWE_TEST, pubKS128) {
    int64_t N_in = 45;
    int64_t N_out = 32;
    int64_t nblimbs_in = 5;
    int64_t nblimbs_out = 6;
    int64_t alpha_bits = nblimbs_in * BITS_PER_LIMBS - 2;
    int64_t limb_prec = limb_precision(alpha_bits);

    BigTorusParams bt_params_in(nblimbs_in);
    BigTorusParams bt_params_out(nblimbs_out);
    TLweParams tlwe_params_in(N_in, bt_params_in);
    TRLweParams trlwe_params_out(N_out, bt_params_out);

    shared_ptr<TLweKey> key_in = tlwe_keygen(tlwe_params_in);
    shared_ptr<TLweKey> key_out = tlwe_keygen(trlwe_params_out);

    cout << "keygen at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<pubKsKey128> ks_key = ks_keygen128(
            trlwe_params_out, tlwe_params_in,
            *key_in, *key_out, alpha_bits);
    cout << "end keygen at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // take a random plaintext
    BigTorus plaintext(bt_params_in);
    random(plaintext);

    // encrypt it
    TLwe ciphertext(tlwe_params_in);
    native_encrypt(ciphertext, plaintext, *key_in);

    //keyswitch it
    TRLwe res(trlwe_params_out);
    cout << "pubks at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    pubKS128(res, ciphertext, *ks_key, limb_prec);
    cout << "end pubks at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // the target is a constant BTPoly of the plaintext
    BigTorusPolynomial phase(N_out, bt_params_out);
    native_phase(phase, res, *key_out, limb_prec * BITS_PER_LIMBS);

    //constant term
    RR::SetPrecision(limb_prec * BITS_PER_LIMBS);
    ASSERT_LE(log2Diff(to_RR(phase.getAT(0)), to_RR(plaintext)), -int64_t(alpha_bits) + 6);
    //other terms
    for (UINT64 i = 1; i < trlwe_params_out.N; i++) {
        RR::SetPrecision(limb_prec * BITS_PER_LIMBS);
        ASSERT_LE(log2Diff(to_RR(phase.getAT(i)), to_RR(0)), -int64_t(alpha_bits) + 6);
    }
}

TEST(TRLWE_TEST, pubKS32) {
    int64_t N_in = 4096;
    int64_t N_out = 4096;
    int64_t nblimbs_in = 2;
    int64_t nblimbs_out = 2;
    int64_t alpha_bits = 80;
    int64_t limb_prec = limb_precision(alpha_bits);

    BigTorusParams bt_params_in(nblimbs_in);
    BigTorusParams bt_params_out(nblimbs_out);
    TLweParams tlwe_params_in(N_in, bt_params_in);
    TRLweParams trlwe_params_out(N_out, bt_params_out);

    shared_ptr<TLweKey> key_in = tlwe_keygen(tlwe_params_in);
    shared_ptr<TLweKey> key_out = tlwe_keygen(trlwe_params_out);

    cout << "keygen at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<pubKsKey32> ks_key = ks_keygen32(
            trlwe_params_out, tlwe_params_in,
            *key_in, *key_out, alpha_bits);
    cout << "end keygen at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // take a random plaintext
    BigTorus plaintext(bt_params_in);
    random(plaintext);

    // encrypt it
    TLwe ciphertext(tlwe_params_in);
    native_encrypt(ciphertext, plaintext, *key_in);

    //keyswitch it
    TRLwe res(trlwe_params_out);
    cout << "pubks at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    pubKS32(res, ciphertext, *ks_key, limb_prec);
    cout << "end pubks at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // the target is a constant BTPoly of the plaintext
    BigTorusPolynomial phase(N_out, bt_params_out);
    native_phase(phase, res, *key_out, limb_prec * BITS_PER_LIMBS);

    //constant term
    RR::SetPrecision(limb_prec * BITS_PER_LIMBS);
    ASSERT_LE(log2Diff(to_RR(phase.getAT(0)), to_RR(plaintext)), -int64_t(alpha_bits) + 2);
    //other terms
    for (UINT64 i = 1; i < trlwe_params_out.N; i++) {
        RR::SetPrecision(limb_prec * BITS_PER_LIMBS);
        ASSERT_LE(log2Diff(to_RR(phase.getAT(i)), to_RR(0)), -int64_t(alpha_bits) + 3);
    }
}
