#include <include/gtest/gtest.h>
#include "../TRLwe.h"
#include "../TRGSW.h"

using namespace std;

TEST(TRGSW_TEST, trgsw_external_product) {
    ASSERT_FALSE(2);


    int64_t N = 64;
    int64_t nblimbs_in = 5;


    int64_t alpha_bits = nblimbs_in * BITS_PER_LIMBS - 2;
    UINT64 out_alpha_bits = alpha_bits;


    BigTorusParams bt_params_in(nblimbs_in);

    TLweParams tlwe_params_in(N, bt_params_in);
    shared_ptr<TLweKey> key = tlwe_keygen(tlwe_params_in);

    TRGSWParams trgswParams(N, bt_params_in);
    TRGSW a(trgswParams);


    TRLweParams trlweParams(N, bt_params_in);
    TRLwe b(trlweParams);



}
