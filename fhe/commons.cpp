#include "commons.h"
#include <NTL/ZZ_limbs.h>

NTL_CLIENT;

void to_fixP(BigFixPRef reps, const NTL::RR &a) {
    const int64_t n = reps.params->torus_params.torus_limbs;
    RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    BigTorusRef ta(reps.limbs_raw, &reps.params->torus_params);
    to_torus(ta, a * pow(2, -(reps.params->plaintext_expo + reps.params->level_expo)));
}

void to_torus(BigTorusRef reps, const NTL::RR &a) {
    assert(abs(a) <= 0.5);
    const int64_t n = reps.params->torus_limbs;
    RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    ZZ az = RoundToZZ((a + 2) * pow(2, n * BITS_PER_LIMBS));
    //cout << "az:" << az << endl;
    const uint64_t *limbs = ZZ_limbs_get(az);
    for (int64_t i = 0; i < n; i++) {
        mpn_copyi(reps.limbs_raw, limbs, n);
    }
}

NTL::RR to_RR(const BigFixPRef &a) {
    BigTorusRef ta(a.limbs_raw, &a.params->torus_params);
    RR reps = to_RR(ta);
    reps *= NTL::pow(to_RR(2), to_RR(a.params->level_expo + a.params->plaintext_expo));
    return reps;
}

NTL::RR to_RR(const BigTorusRef &a) {
    const long n = a.params->torus_limbs;
    NTL::RR pow2m32 = NTL::pow(to_RR(2), to_RR(-32));
    static const uint64_t mask32 = 0xFFFFFFFFul;
    NTL::RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    NTL::RR reps;
    reps = 0;
    for (long i = 0; i < n; i++) {
        uint64_t limb = a.limbs_raw[i];
        reps += (limb & mask32);
        reps *= pow2m32;
        reps += (limb >> 32u);
        reps *= pow2m32;
    }
    if (reps >= 0.5) reps -= 1;
    return reps;
}
