#include <include/gtest/gtest.h>
#include <NTL/RR.h>
#include "../TRLwe.h"
#include "../TRGSW.h"
#include "../BigFFT.h"
#include "../arithmetic.h"

using namespace std;
using namespace NTL;

int64_t *debug_plaintext;
TLweKey *debug_key;

TEST(TRGSW_TEST, trgsw_external_product) {


    int64_t N = 4096;
    //int64_t n = 2 * N;
    int64_t nblimbs_in = 3;


    int64_t alpha_bits = 120; //signed
    int64_t bits_a = 15;
    UINT64 out_alpha_bits = alpha_bits - (32 + int(log2(N)));

    BigTorusParams bt_params_in(nblimbs_in);
    TRGSWParams trgswParams(N, bt_params_in);

    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);
    debug_key = key.get();
    //for (int64_t i = 0; i < N; i++) {
    //    key->key[i] = 0;
    //}

    TRGSW a(trgswParams);

    TRLwe b(trgswParams);

    TRLwe reps(trgswParams);


    int64_t *plaintext_a = new int64_t[a.params.N];
    debug_plaintext = plaintext_a;

    for (UINT64 i = 0; i < a.params.N; i++) {
        //plaintext_a[i] = ((i == 0) ? 1 : 0);
        plaintext_a[i] = (rand() % (1l << bits_a)) - (1l << (bits_a - 1));
    }

    intPoly_encrypt(a, plaintext_a, *key, alpha_bits);


    //verify the encryption of a
    BigTorusPolynomial phaseAC(N, bt_params_in);
    BigTorusPolynomial phaseRef(N, bt_params_in);

    //const BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs_in);

    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < int(a.ell); i++) {
            // the phase of a[j][i]
            native_phase_FFT(phaseAC, a.a[j][i][0], a.a[j][i][1], *key);
            // plaintext * phase(h_j,i)
            zero(phaseRef);
            phaseRef.getAT(0).limbs_end[-i / 2 - 1] = ((i % 2 == 0) ? (1ul << 32ul) : 1ul);
            if (j == 0) {
                mpn_neg(phaseRef.getAT(0).limbs_end - nblimbs_in,
                        phaseRef.getAT(0).limbs_end - nblimbs_in, nblimbs_in);
                fft_external_product(phaseRef, key->key, phaseRef, 1, nblimbs_in);
            }
            fft_external_product(phaseRef, plaintext_a, phaseRef, 16, nblimbs_in);
            //cout << "phase of a[" << j << "][" << i << "]: " << endl;
            for (int k = 0; k < N; k++) {
                //cout << j << "," << i << "," << k << ": "
                //     << log2Diff(phaseRef.getAT(k), phaseAC.getAT(k)) << endl;
                //cout << to_RR(phaseRef.getAT(k)) << endl;
                //cout << to_RR(phaseAC.getAT(k)) << endl;
                ASSERT_LE(log2Diff(phaseRef.getAT(k), phaseAC.getAT(k)), -alpha_bits + 2);
            }
        }
    }

    BigTorusPolynomial plaintext_b(N, bt_params_in);
    random(plaintext_b, nblimbs_in);

    native_encrypt(b, plaintext_b, *key, alpha_bits);

    cout << "external 1 start: " << clock() / double(CLOCKS_PER_SEC) << endl;
    external_product(reps, a, b, out_alpha_bits);
    cout << "external 2 start: " << clock() / double(CLOCKS_PER_SEC) << endl;
    external_product(reps, a, b, out_alpha_bits);
    cout << "external 2 end: " << clock() / double(CLOCKS_PER_SEC) << endl;

    BigTorusPolynomial phase(N, bt_params_in);

    native_phase(phase, reps, *key, alpha_bits);

    BigTorusPolynomial phase2(N, bt_params_in);
    fft_external_product(phase2, plaintext_a, plaintext_b, bits_a, limb_precision(out_alpha_bits));

    for (int64_t i = 0; i < N; i++) {
        ASSERT_LE(log2Diff(phase.getAT(i), phase2.getAT(i)), -out_alpha_bits);
    }
}

/*
TEST(TRGSW_TEST, trgsw_external_product_trivial) {
    int64_t N = 64;
    int64_t nblimbs = 6;
    int64_t ell;

    BigTorusParams bt_params(nblimbs);
    TLweParams tlwe_params_in(N, bt_params);
    TRGSWParams trgswParams(N, bt_params);

    TRLwe b(trgswParams);
    TRGSW a(trgswParams);

    //make a trivial ciphertext in b
    zero(b.a[0]);
    random(b.a[1], nblimbs);

    int64_t plaintextA = (rand() & 0xffff) - 32768;

    //make a trivial ciphertext in a
    a.fft_nlimbs = nblimbs;
    a.bits_a = 16;
    a.ell = ell;

    for (int j=0; j<2; j++) {
        for (int i=0; i<ell; i++) {
            for
        }
    }



}
*/
