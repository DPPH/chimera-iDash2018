#include <iostream>
#include <NTL/LLL.h>
#include "BigFixP.h"
#include <NTL/ZZ_limbs.h>

NTL_CLIENT;

RR to_RR(const BigTorusRef& a) {
    const long n = a.params->torus_limbs;
    RR pow2m32 = NTL::pow(to_RR(2),to_RR(-32));
    static const uint64_t mask32 = 0xFFFFFFFFul;
    RR::SetPrecision(n*BITS_PER_LIMBS+2);
    RR reps; reps=0;
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

RR to_RR(const BigFixPRef& a) {
    BigTorusRef ta(a.limbs_raw, &a.params->torus_params);
    RR reps = to_RR(ta);
    reps *= NTL::pow(to_RR(2), to_RR(a.params->level_expo+a.params->plaintext_expo));
    return reps;
}

void to_torus(BigTorusRef reps, const RR& a) {
    assert(abs(a)<=0.5);
    const int64_t n = reps.params->torus_limbs;
    RR::SetPrecision(n*BITS_PER_LIMBS+2);
    ZZ az = RoundToZZ((a + 2) * pow(2, n*BITS_PER_LIMBS));
    //cout << "az:" << az << endl;
    const uint64_t* limbs = ZZ_limbs_get(az);
    for (int64_t i = 0; i<n; i++) {
        mpn_copyi(reps.limbs_raw, limbs, n);
    }
}

void to_fixP(BigFixPRef reps, const RR& a) {
    const int64_t n = reps.params->torus_params.torus_limbs;
    RR::SetPrecision(n*BITS_PER_LIMBS+2);
    BigTorusRef ta(reps.limbs_raw, &reps.params->torus_params);
    to_torus(ta, a * pow(2, -(reps.params->plaintext_expo+reps.params->level_expo)));
}

int main() {
    BigFixPParams params;
    params.torus_params.torus_limbs = 10;
    params.plaintext_expo = 10;
    params.level_expo = 50;
    BigFixPParams paramsb;
    paramsb.torus_params.torus_limbs = 10;
    paramsb.plaintext_expo = 10;
    paramsb.level_expo = 50;
    BigFixPParams paramsc;
    paramsc.torus_params.torus_limbs = 10;
    paramsc.plaintext_expo = 10;
    paramsc.level_expo = 49;

    //torus load and store to RR
    for (int trial = 0; trial < 100; trial++) {
        RR test1 = random_RR() - 0.5;
        BigTorus testTorus(&params.torus_params);
        to_torus(testTorus, test1);
        RR test2 = to_RR(testTorus);

        RR::SetOutputPrecision(long(640 * log(2) / log(10)));
        //cout << test1 << endl;
        //cout << test2 << endl;
        assert(abs(test1-test2)<=1e-100);
    }
    //fixP load and store to RR
    for (int trial = 0; trial < 100; trial++) {
        RR test1 = 1000 * (random_RR() - 0.5);
        BigFixP testFixp(&params);
        to_fixP(testFixp, test1);
        RR test2 = to_RR(testFixp);

        RR::SetOutputPrecision(long(640 * log(2) / log(10)));
        //cout << test1 << endl;
        //cout << test2 << endl;
        assert(abs(test1-test2)<=1e-100);
    }
    //fixP addition
    for (int trial = 0; trial < 100; trial++) {
        RR testa = 1000 * (random_RR() - 0.5);
        RR testb = 1000 * (random_RR() - 0.5);
        BigFixP testFixpA(&params);
        BigFixP testFixpB(&paramsb);
        BigFixP testFixpC(&paramsc);
        to_fixP(testFixpA, testa);
        to_fixP(testFixpB, testb);
        add(testFixpC, testFixpA, testFixpB, NA);
        RR testc = to_RR(testFixpC);

        RR::SetOutputPrecision(long(640 * log(2) / log(10)));
        cout << testa << endl;
        cout << testb << endl;
        cout << testc << endl;
        assert(abs(testa+testb-testc)<=1e-100);
    }

    return 0;
}
