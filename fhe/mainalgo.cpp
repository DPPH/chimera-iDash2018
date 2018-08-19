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
