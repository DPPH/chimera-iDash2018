#include "mainalgo.h"
#include "BigFFT.h"

NTL_CLIENT;

std::shared_ptr<TRGSWMatrix>
encrypt_S(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits) {
    const int64_t Ns2 = N / 2;
    const int64_t n = N * 2;
    const int64_t rows = plaintext.NumRows();
    const int64_t cols = ceil(plaintext.NumCols() / double(Ns2));

    const int64_t enc_limbs = limb_precision(plaintext_precision_bits + log2(N));
    BigComplex *cslots = new_BigComplex_array(Ns2, enc_limbs);
    BigReal *rcoefs = new_BigReal_array(N, enc_limbs);
    int64_t *icoefs = new int64_t[N];
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, enc_limbs);


    int64_t reps_limbs = limb_precision(alpha_bits);
    int64_t reps_plaintext_exponent = plaintext_precision_bits; //TODO
    int64_t reps_level_exponent = 0; //TODO
    shared_ptr<BigTorusParams> bt_params = make_shared<BigTorusParams>(reps_limbs, reps_plaintext_exponent,
                                                                       reps_level_exponent);
    store_forever(bt_params);
    shared_ptr<TRGSWParams> trgsw_params = make_shared<TRGSWParams>(N, *bt_params);
    store_forever(trgsw_params);
    TRGSWMatrix *reps = new TRGSWMatrix(rows, cols, *trgsw_params);

    for (int64_t r = 0; r < rows; r++) {
        for (int64_t c = 0; c < cols; c++) {
            RR::SetPrecision(enc_limbs * BITS_PER_LIMBS);
            //set the slots at position (r,c)
            for (int64_t slotid = 0; slotid < Ns2; slotid++) {
                int64_t pid = c * Ns2 + slotid;
                if (pid < cols) {
                    to_BigReal(cslots[slotid].real, plaintext[r][pid] / power2_RR(1)); //TODO?
                    zero(cslots[slotid].imag);
                } else {
                    zero(cslots[slotid].real);
                    zero(cslots[slotid].imag);
                }
            }
            //fft
            FFT(rcoefs, cslots, n, powombar);
            //rescale and round
            for (int64_t coefid = 0; coefid < N; coefid++) {
                icoefs[coefid] =
                        to_long(RoundToZZ(to_RR(rcoefs[coefid]) * power2_RR(1 + plaintext_precision_bits))); //TODO
            }
            //encrypt the coefs
            intPoly_encrypt(reps->data[r][c], icoefs, key, alpha_bits);
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
    int64_t actual_plaintext_exponent =  default_plaintext_exponent;
    int64_t actual_level_exponent =      default_level_exponent;
    // If there is a plaintext exponent override, correct the level info
    if (override_plaintext_exponent != int64_t(NA)) {
        int64_t difference = override_plaintext_exponent - actual_plaintext_exponent;
        actual_plaintext_exponent += difference;
        actual_level_exponent -= difference;
    }
    assert_dramatically(actual_level_exponent>0, "impossible to perform the requested multiplication, level too low");
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
    TRLwe shifted_input(shifted_input_trlwe_params);

    //step 2: Do the external product
    int64_t accum_level = target_level_expo;
    int64_t accum_plaintext_expo = actual_plaintext_exponent;
    int64_t accum_alpha_bits = target_level_expo + plaintext_precision_bits;
    int64_t accum_limbs = limb_precision(accum_alpha_bits + 5);
    BigTorusParams accum_bt_params(accum_limbs, accum_plaintext_expo, accum_level);
    TRLweParams accum_trlwe_params(N, accum_bt_params);
    TRLwe accum(accum_trlwe_params);

    int64_t trgsw_alpha_bits = accum_alpha_bits + A.data[0][0].bits_a + log2(N);


    return shared_ptr<TRLWEVector>();
}
