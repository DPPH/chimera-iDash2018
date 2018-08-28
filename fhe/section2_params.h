#ifndef FHE_SECTION2_PARAMS_H
#define FHE_SECTION2_PARAMS_H


#include <cstdint>
#include <string>

namespace section2_params {
    static const int64_t algo_n = 253;
    static const int64_t algo_k = 3;
    static const int64_t algo_m = 10000;
    //
    static const int64_t n_lvl0 = 680;
    static const std::string lvl0_key_filename = "secret_keyset.bin";
    //
    static const int64_t N = 4096;
    static const std::string section2_key_filename = "section2_secret.key";


};


#endif //FHE_SECTION2_PARAMS_H
