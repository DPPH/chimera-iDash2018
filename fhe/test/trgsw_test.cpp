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


TEST(TRGSW_BLINDROTATE_TEST, trgsw_blind_rotate) {
    int64_t N = 4096;
    int64_t n_in = 500;
    int64_t nblimbs = 2;
    int64_t alpha_bits = 120; //signed
    int64_t out_alpha_bits = alpha_bits - (32 + int(log2(N)));

    BigTorusParams bt_params(nblimbs);

    TRGSWParams trgswParams(N, bt_params);

    BigTorusPolynomial phase(N, bt_params);

    TRLwe reps(trgswParams);
    TRGSW *c = new_TRGSW_array(n_in, trgswParams);
    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);
    int64_t b = (rand() % (2 * N));


    int64_t *s = new int64_t[n_in];
    int64_t *a = new int64_t[n_in];
    int64_t power = -b;

    cout << "start encrypt: " << clock() / double(CLOCKS_PER_SEC) << endl;
    for (int i = 0; i < n_in; i++) {
        s[i] = rand() % 2;
        a[i] = rand() % (2 * N);
        binary_encrypt(c[i], s[i], *key, alpha_bits);
        power += s[i] * a[i];
    }
    cout << "end encrypt: " << clock() / double(CLOCKS_PER_SEC) << endl;



    //make a trivial ciphertext in reps
    zero(reps.a[0]);
    random(reps.a[1], nblimbs);
    TRLwe phase_ref(trgswParams);
    copy(phase_ref, reps);


    cout << "start blind rotate: " << clock() / double(CLOCKS_PER_SEC) << endl;
    blind_rotate(reps, b, a, c, n_in, out_alpha_bits);
    cout << "end blind rotate: " << clock() / double(CLOCKS_PER_SEC) << endl;

    native_phase(phase, reps, *key, alpha_bits);
    rotate(phase_ref, phase_ref, power);


    for (int64_t i = 0; i < N; i++) {
        EXPECT_LE(log2Diff(phase.getAT(i), phase_ref.a[1].getAT(i)), -out_alpha_bits);
    }
    delete[] a;
    delete[] s;
    delete_TRGSW_array(n_in, c);
}

