#include <cassert>
#include <omp.h>
#include <fstream>
#include "mainalgo.h"
#include "BigFFT.h"
#include "TRLwe.h"
#include <fstream>
#include <iostream>
#include <cassert>

NTL_CLIENT;

void read_tlwe_key(const char *const filename, int64_t *const key, const int64_t n) {
    static const int32_t LWE_KEY_TYPE_UID = 43;

    ifstream f(filename, ifstream::binary);
    if (f) {

        int32_t type_uid;
        istream_read_binary(f, &type_uid, sizeof(int32_t));
        assert_dramatically(type_uid == LWE_KEY_TYPE_UID);

        int *key_tmp = new int[n];
        istream_read_binary(f, key_tmp, sizeof(int) * n);

        for (int i = 0; i < n; ++i)
            key[i] = (int64_t) key_tmp[i];
        delete[] key_tmp;

        f.close();
    } else {
        fprintf(stderr, "Function %s: cannot open file %s\n", __FUNCTION__, filename);
        srand(42);
        for (int i = 0; i < n; ++i)
            key[i] = rand() % 2;
    }

    printf("Read key:\n");
    for (int i = 0; i < n; ++i)
        printf("%ld ", long(key[i]));
    printf("\n");
}

void read_tlwe_sample(istream &f, int64_t *const ab, const int64_t n, const int64_t modulus) {
    static const int32_t LWE_SAMPLE_TYPE_UID = 42;
    assert_dramatically((modulus & (modulus - 1)) == 0, "modulus must be a power of 2");
    const UINT64 modulus_bits = __builtin_popcount(modulus - 1); //log_2 of modulus
    const UINT64 t64ToModuloShift = 64ul - modulus_bits;

    int32_t type_uid;
    istream_read_binary(f, &type_uid, sizeof(int32_t));
    assert_dramatically(type_uid == LWE_SAMPLE_TYPE_UID);

    istream_read_binary(f, ab, sizeof(int64_t) * (n + 1));

    double variance;
    istream_read_binary(f, &variance, sizeof(double));

    //modulus-rescale from torus64 to modulus
    for (int64_t i = 0; i < n + 1; i++) {
        ab[i] = (UINT64(ab[i]) >> t64ToModuloShift);
        assert(ab[i] >= 0 && ab[i] < modulus);
    }
}

void read_tlwe_samples(const char *const filename, int64_t **const samples, const int64_t nb_samples,
                       const int64_t nb_coefs, const int64_t modulus) {
    ifstream f(filename, ifstream::binary);

    if (f) {
        for (int i = 0; i < nb_samples; ++i)
            read_tlwe_sample(f, samples[i], nb_coefs, modulus);
    } else {
        fprintf(stderr, "Function %s: cannot open file %s\n", __FUNCTION__, filename);
        fprintf(stderr, "WARNING: I'll take random instead\n");
        for (int i = 0; i < nb_samples; ++i) {
            for (int j = 0; j <= nb_coefs; ++j) {
                samples[i][j] = rand() % modulus;
            }
        }
    }

    f.close();
}


std::shared_ptr<TRGSWMatrix>
encrypt_S(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits) {
    const int64_t Ns2 = N / 2;
    const int64_t n = N * 2;
    const int64_t rows = plaintext.NumRows();
    const int64_t cols = ceil(plaintext.NumCols() / double(Ns2));

    //find plaintext maximum
    NTL::RR maxiNorm2Sq;
    maxiNorm2Sq = -1;
    NTL::RR maxiNormInf;
    maxiNormInf = -1;
    for (int r = 0; r < rows; r++) {
        NTL::RR norm2sq;
        for (int c = 0; c < cols; c++) {
            norm2sq = 0;
            for (int k = 0; k < Ns2; k++) {
                if (c * Ns2 + k < plaintext.NumCols()) {
                    RR slot = plaintext[r][c * Ns2 + k];
                    norm2sq += slot * slot;
                    if (abs(slot) > maxiNormInf) maxiNormInf = abs(slot);
                }
            }
            if (norm2sq > maxiNorm2Sq) maxiNorm2Sq = norm2sq;
        }
    }
    //this is the actual plaintext exponent
    long real_plaintext_exponent_a = to_long(ceil(log(maxiNorm2Sq) / 2 / log(2.))); //todo override?
    long real_plaintext_exponent_b = to_long(ceil(log(maxiNormInf) / log(2.)));     //todo override?
    cout << "deduced plaintext exponents: " << real_plaintext_exponent_a << " or " << real_plaintext_exponent_b << endl;
    long real_plaintext_exponent = std::min(real_plaintext_exponent_b, real_plaintext_exponent_a);
    long trgsw_bits_a = plaintext_precision_bits;
    //long trgsw_shift = real_plaintext_exponent - trgsw_bits_a;
    //info: the trgsw will be shifted by -trgsw_shift upon encryption
    //info: upon external product, * the level will decrease by tgrsw_bits_a
    //                             * the plaintext expo increases by real_plaintext_exponent

    //In the first section, we need to compute the FFT of the slots, shift by -trgsw_shift, and round to integers
    // (in practice, we will actually
    // * pre-shift by -plaintext_exponent
    // * do the FFT in the [-1,1] domain
    // * post shift by +plaintext_expo-trgsw_shift=trgsw_bits_a
    const int64_t plain_limbs = limb_precision(plaintext_precision_bits + log2(N));
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, plain_limbs);

    //Then, for the TRGSW encryption itself:
    //  * alpha is the requested error level
    const int64_t reps_limbs = limb_precision(alpha_bits);
    const int64_t reps_plaintext_exponent = 0; //these two are unused (we leave 0: no shift, for test purposes)
    const int64_t reps_level_exponent = 0;     //these two are unused (we leave 0: no shift, for test purposes)
    shared_ptr<BigTorusParams> bt_params = make_shared<BigTorusParams>(reps_limbs, reps_plaintext_exponent,
                                                                       reps_level_exponent);
    store_forever(bt_params);
    shared_ptr<TRGSWParams> trgsw_params = make_shared<TRGSWParams>(N, *bt_params);
    store_forever(trgsw_params);
    TRGSWMatrix *reps = new TRGSWMatrix(rows, cols, *trgsw_params);

#pragma omp parallel for
    for (int64_t r = 0; r < rows; r++) {
        BigComplex *cslots = new_BigComplex_array(Ns2, plain_limbs);
        BigReal *rcoefs = new_BigReal_array(N, plain_limbs);
        int64_t *icoefs = new int64_t[N];

        for (int64_t c = 0; c < cols; c++) {
            RR::SetPrecision(plain_limbs * BITS_PER_LIMBS);
            //set the slots at position (r,c)
            for (int64_t slotid = 0; slotid < Ns2; slotid++) {
                int64_t pid = c * Ns2 + slotid;
                //pre-scale the input by -real_plaintext_exponent
                if (pid < plaintext.NumCols()) {
                    to_BigReal(cslots[slotid].real, plaintext[r][pid] / power2_RR(real_plaintext_exponent));
                    zero(cslots[slotid].imag);
                } else {
                    zero(cslots[slotid].real);
                    zero(cslots[slotid].imag);
                }
            }
            //fft
            FFT(rcoefs, cslots, n, powombar);

            RR::SetPrecision(plain_limbs * BITS_PER_LIMBS);
            //rescale and round
            for (int64_t coefid = 0; coefid < N; coefid++) {
                icoefs[coefid] =
                        to_long(RoundToZZ(to_RR(rcoefs[coefid]) * power2_RR(trgsw_bits_a)));
                //cerr << icoefs[coefid] << " ";
            }
            //cerr << endl;
            //encrypt the coefs
            intPoly_encrypt(reps->data[r][c], icoefs, key, alpha_bits);
            //finally, set the real plaintext exponent and bits_a
            reps->data[r][c].bits_a = trgsw_bits_a;
            reps->data[r][c].plaintext_exponent = real_plaintext_exponent;
        }
        //free resources
        delete_BigComplex_array(Ns2, cslots);
        delete_BigReal_array(N, rcoefs);
        delete[] icoefs;
    }
    //free resources
    return std::shared_ptr<TRGSWMatrix>(reps);
}

std::shared_ptr<TRGSWMatrix>
encrypt_X(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits) {
    const int64_t rows = plaintext.NumRows();
    const int64_t cols = plaintext.NumCols();

    //find plaintext maximum
    NTL::RR maxiNormInf;
    maxiNormInf = -1;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            RR slot = plaintext[r][c];
            if (abs(slot) > maxiNormInf) maxiNormInf = abs(slot);
        }
    }
    //this is the actual plaintext exponent
    long real_plaintext_exponent = to_long(ceil(log(maxiNormInf) / log(2.)));     //todo override?
    cout << "plaintext exponents: " << real_plaintext_exponent << endl;
    long trgsw_bits_a = plaintext_precision_bits;
    const int64_t plain_limbs = limb_precision(plaintext_precision_bits + log2(N));


    //Then, for the TRGSW encryption itself:
    //  * alpha is the requested error level
    int64_t reps_limbs = limb_precision(alpha_bits);
    int64_t reps_plaintext_exponent = 0; //these two are unused (we leave 0: no shift, for test purposes)
    int64_t reps_level_exponent = 0;     //these two are unused (we leave 0: no shift, for test purposes)
    shared_ptr<BigTorusParams> bt_params = make_shared<BigTorusParams>(reps_limbs, reps_plaintext_exponent,
                                                                       reps_level_exponent);
    store_forever(bt_params);
    shared_ptr<TRGSWParams> trgsw_params = make_shared<TRGSWParams>(N, *bt_params);
    store_forever(trgsw_params);
    TRGSWMatrix *reps = new TRGSWMatrix(rows, cols, *trgsw_params);

#pragma omp parallel for
    for (int64_t r = 0; r < rows; r++) {
        for (int64_t c = 0; c < cols; c++) {
            RR::SetPrecision(plain_limbs * BITS_PER_LIMBS);
            int64_t message = to_long(RoundToZZ(plaintext[r][c] * power2_RR(trgsw_bits_a - real_plaintext_exponent)));
            int_encrypt(reps->data[r][c], message, key, alpha_bits);
            //finally, set the real plaintext exponent and bits_a
            reps->data[r][c].bits_a = trgsw_bits_a;
            reps->data[r][c].plaintext_exponent = real_plaintext_exponent;
        }
    }
    return std::shared_ptr<TRGSWMatrix>(reps);
}

std::shared_ptr<TRLWEVector>
encrypt_individual_trlwe(NTL::vec_RR plaintext, const TLweKey &key, int64_t N, int64_t level_expo,
                         int64_t plaintext_expo,
                         int64_t plaintext_precision_bits) {
    const int64_t nbelems = plaintext.length();
    //auto-deduce plaintext exponent
    if (plaintext_expo == int64_t(NA)) {
        //plaintext expo is the smallest exponent s.t. plaintext / 2^plaintext_expo is in [-1,1]
        RR maxPlaintextAbs;
        maxPlaintextAbs = 0;
        for (int i = 0; i < nbelems; i++) {
            if (abs(plaintext[i]) > maxPlaintextAbs) maxPlaintextAbs = abs(plaintext[i]);
        }
        if (maxPlaintextAbs == 0)
            plaintext_expo = -1;
        else
            plaintext_expo = ceil(to_double(log(maxPlaintextAbs)) / log(2.));
        cerr << "in encrypt: deduced plaintext_expo set to " << plaintext_expo << endl;
    }
    //generate the TRLWE parameters
    int64_t enc_alpha_bits = level_expo + plaintext_precision_bits + 1;
    int64_t enc_limbs = limb_precision(enc_alpha_bits);

    shared_ptr<BigTorusParams> enc_bt_params = make_shared<BigTorusParams>(enc_limbs, plaintext_expo, level_expo);
    store_forever(enc_bt_params);
    shared_ptr<TRLweParams> enc_trlwe_params(new TRLweParams(N, *enc_bt_params));
    store_forever(enc_trlwe_params);

    //encrypt each element one by one
    TRLWEVector *reps = new TRLWEVector(nbelems, *enc_trlwe_params);
    for (int i = 0; i < nbelems; i++) {
        fixp_encrypt_number(reps->data[i], plaintext[i], key, plaintext_precision_bits);
    }

    return shared_ptr<TRLWEVector>(reps);
}

std::shared_ptr<TRLWEVector>
vec_mat_prod(const TRLWEVector &v, const TRGSWMatrix &A, int64_t target_level_expo, int64_t override_plaintext_exponent,
             int64_t plaintext_precision_bits) {
    //true plaintext: v * 2^(Lv + pl_v) * A * 2^shift_expo  where plaintext_expo = bits_a + shift_expo.
    //        actual expo:          Lv + pl_v + shift_expo = Lv + pl_v + plaintext_expo - bits_a
    //          new level:          Lv - bits_a
    // ->  new plain expo:          pl_v + pl_A

    //plaintext exponent of the result
    int64_t default_plaintext_exponent = v.data[0].params.fixp_params.plaintext_expo +
                                         A.data[0][0].plaintext_exponent;
    int64_t default_level_exponent = v.data[0].params.fixp_params.level_expo -
                                     A.data[0][0].bits_a;
    int64_t actual_plaintext_exponent = default_plaintext_exponent;
    int64_t actual_level_exponent = default_level_exponent;
    // If there is a plaintext exponent override, correct the level info
    if (override_plaintext_exponent != int64_t(NA)) {
        int64_t difference = override_plaintext_exponent - actual_plaintext_exponent;
        actual_plaintext_exponent += difference;
        actual_level_exponent -= difference;
    }
    assert_dramatically(actual_level_exponent > 0, "impossible to perform the requested multiplication, level too low");
    if (target_level_expo == int64_t(NA)) {
        target_level_expo = actual_level_exponent;
    } else {
        assert_dramatically(target_level_expo > 0,
                            "invalid target level expo requested");
        assert_dramatically(actual_level_exponent >= target_level_expo,
                            "impossible to perform the requested multiplication, level too low");
    }
    //RS to apply to the input ciphertext before calling the external product
    int64_t rsBits = actual_level_exponent - target_level_expo;

    //step 1: RS left shift the input by rsBits
    int64_t N = v.data[0].params.N;
    int64_t shifted_input_level = v.data[0].params.fixp_params.level_expo - rsBits;
    int64_t shifted_input_plaintext_expo = v.data[0].params.fixp_params.plaintext_expo;
    int64_t shifted_input_alpha_bits = shifted_input_level + plaintext_precision_bits;
    int64_t shifted_input_limbs = limb_precision(shifted_input_alpha_bits + 5);
    BigTorusParams shifted_input_bt_params(shifted_input_limbs, shifted_input_plaintext_expo, shifted_input_level);
    TRLweParams shifted_input_trlwe_params(N, shifted_input_bt_params);
    TRLWEVector shifted_input(A.rows, shifted_input_trlwe_params);

    //step 2: Prepare the external product
    int64_t accum_level = target_level_expo;
    int64_t accum_plaintext_expo = actual_plaintext_exponent;
    int64_t accum_alpha_bits = target_level_expo + plaintext_precision_bits;
    int64_t accum_limbs = limb_precision(accum_alpha_bits + 5);
    shared_ptr<BigTorusParams> accum_bt_params = make_shared<BigTorusParams>(accum_limbs, accum_plaintext_expo,
                                                                             accum_level);
    shared_ptr<TRLweParams> accum_trlwe_params = make_shared<TRLweParams>(N, *accum_bt_params);

    int64_t trgsw_alpha_bits = accum_alpha_bits + A.data[0][0].bits_a + log2(N);
    int64_t ell = ceil(trgsw_alpha_bits / double(TRGSWParams::Bgbits));
    assert_dramatically(ell <= int64_t(A.data[0][0].ell));

    store_forever(accum_bt_params);
    store_forever(accum_trlwe_params);
    TRLWEVector *reps = new TRLWEVector(A.cols, *accum_trlwe_params);

    for (int j = 0; j < A.rows; j++) {
        lshift(shifted_input.data[j], v.data[j], rsBits);
    }

    for (int j = 0; j < A.cols; j++) {
        zero(reps->data[j]);
    }
#pragma omp parallel for
    for (int i = 0; i < A.rows; i++) {
        TRLwe accum(*accum_trlwe_params);
        for (int j = 0; j < A.cols; j++) {
            external_product(accum, A.data[i][j], shifted_input.data[i], accum_alpha_bits);
            add(reps->data[j], reps->data[j], accum);
        }
    }

    return shared_ptr<TRLWEVector>(reps);
}

//return a vector with the real part of all slots
NTL::vec_RR decrypt_heaan_packed_trlwe(const TRLWEVector &ciphertext, const TLweKey &key, int64_t length) {
    int64_t Ns2 = ciphertext.data[0].params.N / 2;
    assert(int64_t(ciphertext.length * Ns2) >= length);
    assert(int64_t ((ciphertext.length-1) * Ns2) < length);
    vec_RR reps;
    reps.SetLength(length);
    for (int64_t j = 0; j < ciphertext.length; j++) {
        vec_RR plaintexts = slot_decrypt(ciphertext.data[j], key);
        for (int64_t k = 0; k < Ns2; k++) {
            int64_t idx = j * Ns2 + k;
            if (idx < length) reps[idx] = plaintexts[k];
        }
    }
    return reps;
}

//return a vector with the real part of all slots
NTL::mat_RR decrypt_heaan_packed_trlwe(const TRLweMatrix &ciphertext, const TLweKey &key, int64_t plaintext_cols) {
    int64_t Ns2 = ciphertext.data[0][0].params.N / 2;
    assert(int64_t(ciphertext.cols * Ns2) >= plaintext_cols);
    assert(int64_t ((ciphertext.cols-1) * Ns2) < plaintext_cols);
    mat_RR reps;
    reps.SetDims(ciphertext.rows, plaintext_cols);
    for (int64_t i = 0; i < ciphertext.rows; i++) {
        for (int64_t j = 0; j < ciphertext.cols; j++) {
            vec_RR plaintexts = slot_decrypt(ciphertext.data[i][j], key);
            for (int64_t k = 0; k < Ns2; k++) {
                int64_t idx = j * Ns2 + k;
                if (idx < plaintext_cols) reps[i][idx] = plaintexts[k];
            }
        }
    }
    return reps;
}

NTL::vec_RR decrypt_individual_trlwe(const TRLWEVector &ciphertext, const TLweKey &key, int64_t length) {
    assert(int64_t(ciphertext.length) == length);
    vec_RR reps;
    reps.SetLength(length);
    for (int64_t j = 0; j < length; j++) {
        reps[j] = fixp_decrypt_number(ciphertext.data[j], key);
    }
    return reps;
}

std::shared_ptr<TRLWEVector>
product_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, const TRGSW &rk, int64_t target_level_expo,
                  int64_t override_plaintext_exponent, int64_t plaintext_precision_bits) {

    int64_t N = a.data[0].params.N;
    int64_t length = a.length;
    assert(int64_t(b.length) == length);
    assert(int64_t(b.data[0].params.N) == N);
    assert(int64_t(rk.params.N) == N);

    //tau = tau_a + tau_b
    int64_t tau_a = a.data[0].params.fixp_params.plaintext_expo;
    int64_t tau_b = b.data[0].params.fixp_params.plaintext_expo;
    int64_t tau = tau_a + tau_b;

    //L=min(L_a + L_b) - rho
    int64_t L_a = a.data[0].params.fixp_params.level_expo;
    int64_t L_b = b.data[0].params.fixp_params.level_expo;
    int64_t rho = plaintext_precision_bits;
    int64_t L = std::min(L_a, L_b) - rho;


    if (override_plaintext_exponent != int64_t(NA)) {
        tau = override_plaintext_exponent;
        L = L + tau_a + tau_b - tau;
    }

    int64_t diff = 0;

    if (target_level_expo != int64_t(NA)) {
        assert(target_level_expo < L);
        diff = L - target_level_expo;
        L = target_level_expo;
    }

    int64_t shift_a = diff + L_a - std::min(L_a, L_b);
    int64_t shift_b = diff + L_b - std::min(L_a, L_b);
    int64_t newL_a_b = std::min(L_a, L_b) - diff;

    int64_t new_limbs_a_b = limb_precision(newL_a_b + rho);


    BigTorusParams bt_params_a(new_limbs_a_b, tau_a, newL_a_b);
    BigTorusParams bt_params_b(new_limbs_a_b, tau_b, newL_a_b);
    TRGSWParams params_a(N, bt_params_a);
    TRLwe anew(params_a);
    TRGSWParams params_b(N, bt_params_b);
    TRLwe bnew(params_b);


    int64_t reps_limbs = limb_precision(L + rho);

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau, L);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLWEVector *reps = new TRLWEVector(length, *reps_trlwe_params);

    for (int64_t i = 0; i < length; i++) {
        lshift(anew, a.data[i], shift_a);
        lshift(bnew, b.data[i], shift_b);
        fixp_internal_product(reps->data[i], anew, bnew, rk, rho);
    }

    return shared_ptr<TRLWEVector>(reps);
}


std::shared_ptr<TRLWEVector> substract_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, int64_t target_level_expo,
                                                 int64_t override_plaintext_exponent,
                                                 int64_t plaintext_precision_bits) {

    int64_t N = a.data[0].params.N;
    int64_t length = a.length;
    assert(int64_t(b.length) == length);
    assert(int64_t(b.data[0].params.N) == N);
    assert(int64_t(b.data[0].params.N) == N);

    int64_t rho = plaintext_precision_bits;

    //tau=max(tau_a, tau_b) + 1
    const int64_t tau_a = a.data[0].params.fixp_params.plaintext_expo; //50   tau_a=83   L_a= 17
    const int64_t tau_b = b.data[0].params.fixp_params.plaintext_expo; //-50  tau_b=70   L_b= 13
    int64_t tau = std::max(tau_a, tau_b) + 1; //51  tau= 84

    //L=min(L_a + tau_a, L_b+ tau_b) - tau
    int64_t L_a = a.data[0].params.fixp_params.level_expo; //70  17
    int64_t L_b = b.data[0].params.fixp_params.level_expo; //70 13

    if (tau_a - tau_b > rho) { return copy(a); }  //17
    else if (tau_b - tau_a > rho) { return neg(b); }


    int64_t L = std::min(L_a + tau_a, L_b + tau_b) - tau;  //83 -84 = -1
    assert(L >= 0);

    if (override_plaintext_exponent != int64_t(NA)) {
        int64_t difftau = override_plaintext_exponent - tau;
        tau = override_plaintext_exponent; // tau + difftau
        L = L - difftau;
    }


    if (target_level_expo != int64_t(NA)) {
        assert(target_level_expo <= L);
        L = target_level_expo;
    }

    int64_t reps_limbs = limb_precision(L + rho);

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau, L);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLWEVector *reps = new TRLWEVector(length, *reps_trlwe_params);

    assert(int64_t(b.length) == length);
    for (int64_t i = 0; i < length; i++) {

        fixp_sub(reps->data[i], a.data[i], b.data[i], L + rho);
    }

    return shared_ptr<TRLWEVector>(reps);
}


std::shared_ptr<TRLWEVector> add_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, int64_t target_level_expo,
                                           int64_t override_plaintext_exponent,
                                           int64_t plaintext_precision_bits) {

    int64_t N = a.data[0].params.N;
    int64_t length = a.length;
    assert(int64_t(b.length) == length);
    assert(int64_t(b.data[0].params.N) == N);
    assert(int64_t(b.data[0].params.N) == N);

    int64_t rho = plaintext_precision_bits;

    //tau=max(tau_a, tau_b) + 1
    int64_t tau_a = a.data[0].params.fixp_params.plaintext_expo; //50   tau_a=83   L_a= 17
    int64_t tau_b = b.data[0].params.fixp_params.plaintext_expo; //-50  tau_b=70   L_b= 13
    int64_t tau = std::max(tau_a, tau_b) + 1; //51  tau= 84

    //L=min(L_a + tau_a, L_b+ tau_b) - tau
    int64_t L_a = a.data[0].params.fixp_params.level_expo; //70  17
    int64_t L_b = b.data[0].params.fixp_params.level_expo; //70 13

    if (tau_a - tau_b > rho) { return copy(a); }  //17
    else if (tau_b - tau_a > rho) { return neg(b); }


    int64_t L = std::min(L_a + tau_a, L_b + tau_b) - tau;  //83 -84 = -1
    assert(L >= 0);

    if (override_plaintext_exponent != int64_t(NA)) {
        tau = override_plaintext_exponent;
        L = L + tau_a + tau_b - tau;
    }


    if (target_level_expo != int64_t(NA)) {
        assert(target_level_expo < L);
        L = target_level_expo;
    }

    int64_t reps_limbs = limb_precision(L + rho);

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau, L);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLWEVector *reps = new TRLWEVector(length, *reps_trlwe_params);

    assert(int64_t(b.length) == length);
    for (int64_t i = 0; i < length; i++) {

        fixp_add(reps->data[i], a.data[i], b.data[i], L + rho);
    }

    return shared_ptr<TRLWEVector>(reps);
}

std::shared_ptr<TRLWEVector> copy(const TRLWEVector &a) {
    int64_t reps_limbs = a.data[0].params.fixp_params.torus_limbs;
    int64_t tau = a.data[0].params.fixp_params.plaintext_expo;
    int64_t L = a.data[0].params.fixp_params.level_expo;
    int64_t N = a.data[0].params.N;
    int64_t length = a.length;

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau, L);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLWEVector *reps = new TRLWEVector(length, *reps_trlwe_params);

    for (int64_t i = 0; i < length; i++) {

        copy(reps->data[i], a.data[i]);
    }

    return shared_ptr<TRLWEVector>(reps);
}

std::shared_ptr<TRLWEVector> neg(const TRLWEVector &a) {
    int64_t reps_limbs = a.data[0].params.fixp_params.torus_limbs;
    int64_t tau = a.data[0].params.fixp_params.plaintext_expo;
    int64_t L = a.data[0].params.fixp_params.level_expo;
    int64_t N = a.data[0].params.N;
    int64_t length = a.length;

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau, L);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLWEVector *reps = new TRLWEVector(length, *reps_trlwe_params);

    for (int64_t i = 0; i < length; i++) {
        neg(reps->data[i], a.data[i]);
    }

    return shared_ptr<TRLWEVector>(reps);

}

std::shared_ptr<TRLWEVector>
compute_w(const TRLWEVector &p, const TRGSW &rk, int64_t target_level_expo, int64_t override_plaintext_exponent,
          int64_t plaintext_precision_bits) {

    std::shared_ptr<TRLWEVector> pp = product_ind_TRLWE(p, p, rk, NA, NA,
                                                        plaintext_precision_bits);
    std::shared_ptr<TRLWEVector> resp = substract_ind_TRLWE(p, *pp, target_level_expo, override_plaintext_exponent,
                                                            plaintext_precision_bits);
    return shared_ptr<TRLWEVector>(resp);
}

#ifdef DEBUG_COMPUTE_A
extern NTL::mat_RR debug_S;
extern NTL::mat_RR debug_X;
extern NTL::vec_RR debug_W;
extern shared_ptr<TLweKey> debug_key;
#endif

std::shared_ptr<TRLweMatrix>
compute_A(const TRGSWMatrix &X, const TRGSWMatrix &S, const TRLWEVector &W, int64_t target_level_expo,
          int64_t override_plaintext_exponent, int64_t plaintext_precision_bits) {

    int64_t n = W.length;
    assert(X.rows == n);
    assert(S.rows == n);

    int64_t L = W.data[0].params.fixp_params.level_expo - S.data[0][0].bits_a - X.data[0][0].bits_a;
    int64_t tau = W.data[0].params.fixp_params.plaintext_expo + S.data[0][0].plaintext_exponent +
                  X.data[0][0].plaintext_exponent;
    const int64_t N = W.data[0].params.N;

    int64_t L_temp = W.data[0].params.fixp_params.level_expo - S.data[0][0].bits_a;
    int64_t tau_temp = W.data[0].params.fixp_params.plaintext_expo + S.data[0][0].plaintext_exponent;
    int64_t alpha_bits_temp = L_temp + plaintext_precision_bits + 10;
    int64_t temp_limbs = limb_precision(alpha_bits_temp + 10);

    int64_t tau_reps;
    int64_t L_reps;
    if (override_plaintext_exponent != int64_t(NA)) {
        int64_t tau_diff = override_plaintext_exponent - tau;
        tau_reps = override_plaintext_exponent; //tau_temp + tau_diff
        L_reps = L - tau_diff;
    } else {
        tau_reps = tau;
        L_reps = L;
    }
    if (target_level_expo != int64_t(NA)) {
        int64_t L_diff = L_reps - target_level_expo;
        L_reps = target_level_expo;
        assert_dramatically(L_diff == 0, "requested level too high");
        //TODO: allow >=0 and do a lshift
    }
    int64_t reps_alpha_bits = L_reps + plaintext_precision_bits;
    int64_t reps_limbs = limb_precision(reps_alpha_bits);

    std::shared_ptr<BigTorusParams> temp_bt_params = make_shared<BigTorusParams>(temp_limbs, tau_temp, L_temp);
    store_forever(temp_bt_params);

    std::shared_ptr<TRLweParams> temp_trlwe_params = make_shared<TRLweParams>(N, *temp_bt_params);
    store_forever(temp_trlwe_params);

    std::shared_ptr<BigTorusParams> reps_bt_params = make_shared<BigTorusParams>(reps_limbs, tau_reps, L_reps);
    store_forever(reps_bt_params);

    std::shared_ptr<TRLweParams> reps_trlwe_params = make_shared<TRLweParams>(N, *reps_bt_params);
    store_forever(reps_trlwe_params);

    TRLweMatrix *reps = new TRLweMatrix(X.cols, S.cols, *reps_trlwe_params);

    //We'll use one temp trlwe per thread.
    const int NUM_THREADS = omp_get_max_threads();
    cout << "Compute A will use: " << NUM_THREADS << "threads" << endl;
    TRLWEVector temp(NUM_THREADS, *temp_trlwe_params);
    TRLWEVector temp2(NUM_THREADS, *reps_trlwe_params);

    for (int i = 0; i < X.cols; i++) {
        for (int j = 0; j < S.cols; j++) {
            zero(reps->data[i][j]);
        }
    }
#ifdef DEBUG_COMPUTE_A
    BigTorusParams debug_Sdec1_params(S.data[0][0].fft_nlimbs, 0, 0);
    BigTorusPolynomial debug_Sdec1(N, debug_Sdec1_params);
    BigComplex* debug_Sdec2 = new_BigComplex_array(N/2, S.data[0][0].fft_nlimbs);
#endif
    for (int i = 0; i < X.cols; i++) {
        for (int j = 0; j < S.cols; j++) {
#pragma omp parallel for //here, we choose to parallelize the loop in n
            for (int k = 0; k < n; k++) {
                int threadNo = omp_get_thread_num();
                external_product(temp.data[threadNo], S.data[k][j], W.data[k], alpha_bits_temp);
#ifdef DEBUG_COMPUTE_A
                //debug: temp must be equal to W[k] * S[k][...]
                vec_RR debug_tmp = slot_decrypt(temp, *debug_key);
                vec_RR debug_expect; debug_expect.SetLength(N/2);
                //vec_RR debug_Sdec = slot_decrypt(S.data[k][j].a[1][0]);
                native_phase_FFT(debug_Sdec1, S.data[k][j].a[1][0][0], S.data[k][j].a[1][0][1], *debug_key);
                iFFT(debug_Sdec2, debug_Sdec1);

                for (int kk=0; kk<N/2; kk++) {
                    debug_expect[kk] = debug_W[k]*debug_S[k][j*N/2+kk];
                    RR decS = to_RR(debug_Sdec2[kk].real) * power2_RR(S.data[k][j].plaintext_exponent + 32 - S.data[k][j].bits_a);
                    cout << "s..: " << decS << endl << "exs: " << debug_S[k][j*N/2+kk] << endl;
                    cout << "w..: " << slot_decrypt(W.data[k], *debug_key)[0] << endl << "exw: " << debug_W[k] << endl;
                    cout << "dec: " << debug_tmp[kk] << endl << "exp: " << debug_expect[kk] << endl;
                }
#endif
                external_product(temp2.data[threadNo], X.data[k][i], temp.data[threadNo], reps_alpha_bits);
#pragma omp critical
                {
                    fixp_add(reps->data[i][j], reps->data[i][j],
                             temp2.data[threadNo]); //if temp2 and reps are not with the same params we should use fixp_add
                }
            }
        }
    }

    return shared_ptr<TRLweMatrix>(reps);
}

