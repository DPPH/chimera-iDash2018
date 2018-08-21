#include <cassert>
#include "mainalgo.h"
#include "BigFFT.h"
#include "TRLwe.h"

NTL_CLIENT;

std::shared_ptr<TRGSWMatrix>
encrypt_S(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits) {
    const int64_t Ns2 = N / 2;
    const int64_t n = N * 2;
    const int64_t rows = plaintext.NumRows();
    const int64_t cols = ceil(plaintext.NumCols() / double(Ns2));

    //find plaintext maximum
    NTL::RR maxi;
    maxi = -1;
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (abs(plaintext[r][c]) > maxi)
                maxi = abs(plaintext[r][c]);
        }
    }
    long rescale_bits = long(ceil(to_double(log(maxi)) / log(2.)));
    cerr << "rescale_bits: " << rescale_bits << endl;

    const int64_t enc_limbs = limb_precision(plaintext_precision_bits + log2(N));
    BigComplex *cslots = new_BigComplex_array(Ns2, enc_limbs);
    BigReal *rcoefs = new_BigReal_array(N, enc_limbs);
    int64_t *icoefs = new int64_t[N];
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, enc_limbs);

    int64_t reps_limbs = limb_precision(alpha_bits);
    int64_t reps_plaintext_exponent = 0; //these two are unused
    int64_t reps_level_exponent = 0;     //these two are unused
    shared_ptr<BigTorusParams> bt_params = make_shared<BigTorusParams>(reps_limbs, reps_plaintext_exponent,
                                                                       reps_level_exponent);
    store_forever(bt_params);
    shared_ptr<TRGSWParams> trgsw_params = make_shared<TRGSWParams>(N, *bt_params);
    store_forever(trgsw_params);
    TRGSWMatrix *reps = new TRGSWMatrix(rows, cols, *trgsw_params);

    NTL::mat_RR coefs_matrix;
    coefs_matrix.SetDims(rows, cols * N);
    for (int64_t r = 0; r < rows; r++) {
        for (int64_t c = 0; c < cols; c++) {
            RR::SetPrecision(enc_limbs * BITS_PER_LIMBS);
            //set the slots at position (r,c)
            for (int64_t slotid = 0; slotid < Ns2; slotid++) {
                int64_t pid = c * Ns2 + slotid;
                if (pid < plaintext.NumCols()) {
                    to_BigReal(cslots[slotid].real, plaintext[r][pid] / power2_RR(rescale_bits));
                    zero(cslots[slotid].imag);
                } else {
                    zero(cslots[slotid].real);
                    zero(cslots[slotid].imag);
                }
            }
            //fft
            FFT(rcoefs, cslots, n, powombar);
            //store it to the RR coefs matrix
            for (int64_t coefid = 0; coefid < N; coefid++) {
                coefs_matrix[r][c * N + coefid] = to_RR(rcoefs[coefid]) * power2_RR(rescale_bits);
            }
        }
    }

    int64_t trgsw_plaintext_expo = 0;
    int64_t trgsw_bits_a = 0;
    int64_t trgsw_bits_shift = 0;
    //auto-deduce plaintext exponent
    {
        //plaintext expo is the smallest exponent s.t. plaintext / 2^plaintext_expo is in [-1,1]
        RR maxAbsPlaintext;
        maxAbsPlaintext = 0;
        RR maxPlaintextNorm2;
        maxPlaintextNorm2 = 0;
        RR sumPlaintextSq;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sumPlaintextSq = 0;
                for (int k = 0; k < N; k++) {
                    sumPlaintextSq += pow(plaintext[i][j], to_RR(2));
                    if (abs(plaintext[i][j]) > maxAbsPlaintext)
                        maxAbsPlaintext = abs(plaintext[i][j]);
                }
                if (sumPlaintextSq > maxPlaintextNorm2) {
                    maxPlaintextNorm2 = sumPlaintextSq;
                }
            }
        }
        if (maxPlaintextNorm2 == 0) {
            trgsw_plaintext_expo = 0;
            trgsw_bits_a = 0;
        } else {
            //we shift the trgsw by this amount
            trgsw_bits_shift = plaintext_precision_bits - ceil(to_double(log(maxAbsPlaintext)) / log(2.));
            //new norm (level increase)
            trgsw_bits_a = ceil(to_double(log(sumPlaintextSq)) / 2. * log(2.)) + trgsw_bits_shift;
            trgsw_plaintext_expo = -trgsw_bits_shift;
        }
        cerr << "in encrypt: deduced plaintext_expo set to " << trgsw_plaintext_expo << endl;
        cerr << "in encrypt: deduced plaintext_bits_a set to " << trgsw_bits_a << endl;
        cerr << "in encrypt: deduced bits_shift set to " << trgsw_bits_shift << endl;
    }


    for (int64_t r = 0; r < rows; r++) {
        for (int64_t c = 0; c < cols; c++) {
            RR::SetPrecision(enc_limbs * BITS_PER_LIMBS);
            //rescale and round
            for (int64_t coefid = 0; coefid < N; coefid++) {
                icoefs[coefid] =
                        to_long(RoundToZZ(coefs_matrix[r][c * N + coefid] * power2_RR(trgsw_bits_shift)));
                cerr << icoefs[coefid] << " ";
            }
            cerr << endl;
            //encrypt the coefs
            intPoly_encrypt(reps->data[r][c], icoefs, key, alpha_bits);
            reps->data[r][c].bits_a = trgsw_bits_a;
            reps->data[r][c].plaintext_exponent = trgsw_plaintext_expo;
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

    //encrypt each element one by one
    TRLWEVector *reps = new TRLWEVector(nbelems, *enc_trlwe_params);
    for (int i = 0; i < nbelems; i++) {
        fixp_encrypt_number(reps->data[i], plaintext[i], key, plaintext_precision_bits);
    }

    return shared_ptr<TRLWEVector>(reps);
}

std::shared_ptr<TRLWEVector>
mat_vec_prod(const TRLWEVector &v, const TRGSWMatrix &A, int64_t target_level_expo, int64_t override_plaintext_exponent,
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
    TRLWEVector shifted_input(A.cols, shifted_input_trlwe_params);

    //step 2: Prepare the external product
    int64_t accum_level = target_level_expo;
    int64_t accum_plaintext_expo = actual_plaintext_exponent;
    int64_t accum_alpha_bits = target_level_expo + plaintext_precision_bits;
    int64_t accum_limbs = limb_precision(accum_alpha_bits + 5);
    shared_ptr<BigTorusParams> accum_bt_params = make_shared<BigTorusParams>(accum_limbs, accum_plaintext_expo,
                                                                             accum_level);
    shared_ptr<TRLweParams> accum_trlwe_params = make_shared<TRLweParams>(N, *accum_bt_params);
    TRLwe accum(*accum_trlwe_params);

    int64_t trgsw_alpha_bits = accum_alpha_bits + A.data[0][0].bits_a + log2(N);
    int64_t ell = ceil(trgsw_alpha_bits / double(TRGSWParams::Bgbits));
    assert(ell <= int64_t(A.data[0][0].ell));

    store_forever(accum_bt_params);
    store_forever(accum_trlwe_params);
    TRLWEVector *reps = new TRLWEVector(A.cols, *accum_trlwe_params);

    for (int j = 0; j < A.cols; j++) {
        lshift(shifted_input.data[j], v.data[j], rsBits);
    }

    for (int j = 0; j < A.cols; j++) {
        zero(reps->data[j]);
    }

    for (int i = 0; i < A.rows; i++) {
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

    int64_t length = a.length;
    assert(int64_t(b.length) == length);


    //tau=tau_a + tau_b
    int64_t tau_a = a.data[0].params.fixp_params.plaintext_expo;
    int64_t tau_b = b.data[0].params.fixp_params.plaintext_expo;
    int64_t tau = tau_a + tau_b;


    //L=min(L_a + L_b) - rho
    int64_t L_a = a.data[0].params.fixp_params.level_expo;
    int64_t L_b = b.data[0].params.fixp_params.level_expo;
    int64_t rho = plaintext_precision_bits;
    int64_t L = std::min(L_a, L_b) - rho;
    int64_t alpha_bits = L + rho;
    int64_t reps_limbs = limb_precision(alpha_bits);


    BigTorusParams reps_bt_params(reps_limbs, tau, L);
    TRLweParams reps_trlwe_params(length, reps_bt_params);
    TRLWEVector *reps = new TRLWEVector(length, reps_trlwe_params);

    for (int64_t i = 0; i < length; i++) {
        fixp_internal_product(reps->data[i], a.data[i], b.data[i], rk, rho);
    }

    return shared_ptr<TRLWEVector>(reps);
}

std::shared_ptr<TRLWEVector> substract_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, int64_t target_level_expo,
                                                 int64_t override_plaintext_exponent,
                                                 int64_t plaintext_precision_bits) {

    int64_t length = a.length;
    assert(int64_t(b.length) == length);

    int64_t rho = plaintext_precision_bits;

    //tau=max(tau_a, tau_b) + 1
    int64_t tau_a = a.data[0].params.fixp_params.plaintext_expo;
    int64_t tau_b = b.data[0].params.fixp_params.plaintext_expo;
    int64_t tau = std::max(tau_a, tau_b) + 1;

    //L=min(L_a + tau_a, L_b+ tau_b) - tau
    int64_t L_a = a.data[0].params.fixp_params.level_expo;
    int64_t L_b = b.data[0].params.fixp_params.level_expo;
    int64_t L = std::min(L_a + tau_a, L_b + tau_b) - tau;

    int64_t alpha_bits = L + rho;
    int64_t reps_limbs = limb_precision(alpha_bits);


    BigTorusParams reps_bt_params(reps_limbs, tau, L);
    TRLweParams reps_trlwe_params(length, reps_bt_params);
    TRLWEVector *reps = new TRLWEVector(length, reps_trlwe_params);

    for (int64_t i = 0; i < length; i++) {
        sub(reps->data[i], a.data[i], b.data[i]);
    }
    return shared_ptr<TRLWEVector>(reps);
}
