#ifndef FHE_SECTION2_PARAMS_H
#define FHE_SECTION2_PARAMS_H


#include <cstdint>
#include <string>
#include "commons.h"

namespace section2_params {
    static const int64_t algo_n = 245;
    static const int64_t algo_k = 3;
    static const int64_t algo_m = 10000;
    //
    static const int64_t default_plaintext_precision = 16;
    //
    static const int64_t n_lvl0 = 680;    //TODO: synchronize with section 1
    static const std::string lvl0_key_filename = "secret_keyset.bin";
    //
    static const int64_t N = 4096;
    static const std::string section2_key_filename = "section2_secret.key";
    //
    static const std::string section1_2_bk_filename = "bk.key";
    static const int64_t bk_alpha_bits = 115; //signed
    static const int64_t bk_nblimbs = limb_precision(bk_alpha_bits);
    //
    static const std::string section1_2_ks_filename = "ks.key";
    static const int64_t ks_nblimbs_out = 2;
    static const int64_t ks_out_alpha_bits = 80;
    //
    static const std::string section2_rk_filename = "rk.key";
    static const int64_t rk_alpha_bits = 110;
    static const int64_t rk_nblimbs = limb_precision(rk_alpha_bits);
    //
    // input tlwe during mega-bootstrap (mapping N/2 -> 2^test_vector_plaintext_expo)
    static const std::string lvl0_ciphertext_filename = "X_beta.ctxt";
    static const int64_t lvl0_ciphertext_modulus = 2 * N; //lvl0_ciphertext are Torus64 and will be read as modulo 2N
    static const int64_t lvl0_ciphertext_plaintext_expo = 2; //the test vector is between -4 and 4
    //TODO: synchronize with section 1
    //
    static const std::string p_lvl4_filename = "p_lvl4.bin";
    static const int64_t p_level = 64;
    static const int64_t p_plaintext_expo = 0;
    static const int64_t p_alpha_bits = p_level + default_plaintext_precision;
    static const int64_t p_limbs = limb_precision(p_alpha_bits);

    static const std::string w_lvl3_filename = "w_lvl3.bin";
    static const int64_t w_level = 49;
    static const int64_t w_plaintext_expo = -2;
    static const int64_t w_alpha_bits = w_level + default_plaintext_precision;
    static const int64_t w_limbs = limb_precision(w_alpha_bits);


};


#endif //FHE_SECTION2_PARAMS_H
