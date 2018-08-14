#include "TRGSW.h"

TRGSW::TRGSW(const TRGSWParams &params) : params(params) {

    for (UINT64 i = 0; i < params.ell; i++) {
        for (UINT64 j = 0; i < 2; i++) {

            a[i][j] = new_BigComplex_array(params.N / 2, params.fixp_params.torus_limbs); //TODO
        }
    }

}

TRGSW::~TRGSW() {
    for (UINT64 i = 0; i < params.ell; i++) {
        for (UINT64 j = 0; i < 2; i++) {

            delete_BigComplex_array(params.N / 2, a[i][j]);
        }
    }
}


TRGSWParams::TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params) : TRLweParams(N, fixp_params) {

}


void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits) {
    abort();
}

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 alpha_bits) {
    abort();
}
