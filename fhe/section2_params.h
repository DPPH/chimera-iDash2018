#ifndef FHE_SECTION2_PARAMS_H
#define FHE_SECTION2_PARAMS_H


#include <cstdint>
#include <string>
#include "commons.h"

namespace section2_params {
    static const int64_t algo_n = 245;
    static const int64_t algo_k = 3;
    static const int64_t algo_m = 10643;
    //
    static const int64_t default_plaintext_precision = 14;
    //
    static const int64_t n_lvl0 = 50;    //TODO: synchronize with section 1
    static const std::string lvl0_key_filename = "secret_keyset.bin";
    //
    static const int64_t N = 4096; //TODO
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
    static const int64_t p_level = 66;
    static const int64_t p_plaintext_expo = 0;
    static const int64_t p_alpha_bits = p_level + default_plaintext_precision;
    static const int64_t p_limbs = limb_precision(p_alpha_bits);

    static const std::string w_lvl3_filename = "w_lvl3.bin";
    static const int64_t w_level = 54; //80-2*14+2
    static const int64_t w_plaintext_expo = -2;
    static const int64_t w_alpha_bits = w_level + default_plaintext_precision;
    static const int64_t w_limbs = limb_precision(w_alpha_bits);

    static const std::string y_lvl2_filename = "y_lvl2.bin";
    static const int64_t y_level = 38;
    static const int64_t y_plaintext_expo = 0;

    static const std::string S_lvl3_filename = "S_lvl3.bin";
    static const int64_t S_level = 52;
    static const int64_t S_plaintext_expo = 0;
    static const int64_t S_alpha_bits = S_level + 32 + 5;


    static const std::string X_lvl2_filename = "X_lvl2.bin";
    static const int64_t X_level = 38;
    static const int64_t X_plaintext_expo = 0;
    static const int64_t X_alpha_bits = X_level + 32 + 5;

    static const std::string numerator_lvl0_filename = "numerator_lvl0.bin";
    static const int64_t numerator_level = 1;
    static const int64_t numerator_plaintext_expo = 6; //TODO

    static const std::string A_lvl1_filename = "A_lvl1.bin";
    static const int64_t A_level = 21;
    static const int64_t A_plaintext_expo = 2; //TODO


    static const std::string denominator_lvl0_filename = "denominator_lvl0.bin";
    static const int64_t denominator_level = 1;
    static const int64_t denominator_plaintext_expo = 6; //TODO
    static const int64_t denominator_alpha_bits = denominator_level + default_plaintext_precision;
    static const int64_t denominator_limbs = limb_precision(denominator_alpha_bits);

};


#endif //FHE_SECTION2_PARAMS_H
