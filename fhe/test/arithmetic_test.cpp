#include "gtest/gtest.h"
#include <NTL/LLL.h>
#include "../commons.h"
#include "../arithmetic.h"

NTL_CLIENT;

double log2Diff(const RR &a, const RR &b) {
    if (a == b) return -INFINITY;
    return to_double(log(abs(a - b))) / log(2.);
}

RR centerMod1(const RR &a) {
    return a - to_RR(RoundToZZ(a));
}

double log2DiffMod1(const RR &a, const RR &b) {
    RR dist = centerMod1(a - b);
    if (dist == 0) return -INFINITY;
    return to_double(log(abs(dist))) / log(2.);
}

TEST(BIGTORUS_ARITHMETIC, convRR) {
    for (int torus_limbs: {5, 10, 15}) {
        BigTorusParams params(torus_limbs);

        //torus load and store to RR
        for (int trial = 0; trial < 100; trial++) {
            RR test1 = random_RR() - 0.5;
            BigTorus testTorus(params);
            to_torus(testTorus, test1);
            RR test2 = fixp_to_RR(testTorus);

            RR::SetOutputPrecision(long(640 * log(2) / log(10)));
            //cout << test1 << endl;
            //cout << test2 << endl;
            ASSERT_LE(log2Diff(test1, test2), -torus_limbs * BITS_PER_LIMBS);
        }
    }
}

TEST(BIGTORUS_ARITHMETIC, convFixP) {
    for (int64_t torus_limbs : {5, 10}) {
        for (int64_t plaintext_expo : {10, 20}) {
            for (int64_t level_expo : {20, 63, 128}) {
                BigTorusParams params(torus_limbs, plaintext_expo, level_expo);
                int64_t torus_prec = torus_limbs * BITS_PER_LIMBS - plaintext_expo - level_expo;

                //fixP load and store to RR
                for (int trial = 0; trial < 100; trial++) {
                    RR test1 = 1000 * (random_RR() - 0.5);
                    BigTorus testFixp(params);
                    to_fixP(testFixp, test1);
                    RR test2 = fixp_to_RR(testFixp);

                    RR::SetOutputPrecision(long(640 * log(2) / log(10)));
                    //cout << test1 << endl;
                    //cout << test2 << endl;
                    ASSERT_LE(log2Diff(test1, test2), -torus_prec);
                }
            }
        }
    }
}

TEST(BIGTORUS_ARITHMETIC, addFixP) {
    BigTorusParams params(10, 10, 128);
    BigTorusParams paramsb(11, 10, 128);
    BigTorusParams paramsc(10, 10, 64);

    //fixP addition
    for (int trial = 0; trial < 100; trial++) {
        RR testa = 1000 * (random_RR() - 0.5);
        RR testb = 1000 * (random_RR() - 0.5);
        BigTorus testFixpA(params);
        BigTorus testFixpB(paramsb);
        BigTorus testFixpC(paramsc);
        to_fixP(testFixpA, testa);
        to_fixP(testFixpB, testb);
        fixp_add(testFixpC, testFixpA, testFixpB, 640);
        RR testc = fixp_to_RR(testFixpC);

        RR::SetOutputPrecision(long(640 * log(2) / log(10)));
        //cout << testa << endl;
        //cout << testb << endl;
        //cout << testa+testb << endl;
        //cout << testc << endl;
        ASSERT_LE(log2Diff(testa + testb, testc), -500);
    }
}

#include <NTL/ZZ_limbs.h>

TEST(BIGTORUS_ARITHMETIC, submul128) {
    for (int64_t limb_prec:  {1, 3, 5, 7}) {
        for (int64_t out_nblimbs: {limb_prec, limb_prec + 1, limb_prec + 3}) {
            for (int64_t in_nblimbs: {limb_prec + 2, limb_prec + 3, limb_prec + 5}) {
                for (bool positive: {true, false}) {


                    BigTorusParams in_params(in_nblimbs);
                    BigTorusParams out_params(out_nblimbs);

                    BigTorus a(in_params);
                    RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
                    RR aa = (random_RR() - 0.5);
                    to_torus(a, aa);

                    RR::SetOutputPrecision(100);

                    //cout << endl;
                    //cout << "aa: " << aa << endl;

                    BigTorus out(out_params);
                    RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
                    RR oout = (random_RR() - 0.5);
                    to_torus(out, oout);

                    //cout << "oout: " << oout << endl;


                    //test positive first
                    ZZ coef_zz = RandomBits_ZZ(127);
                    __int128 coef = *(__int128 *) ZZ_limbs_get(coef_zz);
                    if (!positive) {
                        coef = -coef;
                        coef_zz = -coef_zz;
                    }

                    //cout << "coef_zz: " << coef_zz << endl;

                    subMulS128(out, coef, a, limb_prec);

                    //test that out = oldout - a * coef modulo 1
                    RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
                    RR target = centerMod1(oout - to_RR(coef_zz) * aa);
                    RR actual = to_RR(out);

                    //cout << "actual: " << actual << endl;
                    //cout << "target: " << target << endl;

                    RR distance = centerMod1(target - actual);
                    //cout << "distance: " << distance << endl;

                    ASSERT_LE(log2DiffMod1(target, actual), -limb_prec * BITS_PER_LIMBS + 1);
                }
            }
        }
    }
}

TEST(BIGTORUS_ARITHMETIC, mulS64) {
    for (int64_t limb_prec:  {1, 3, 5, 7}) {
        for (int64_t out_nblimbs: {limb_prec, limb_prec + 1, limb_prec + 3}) {
            int64_t in_nblimbs = out_nblimbs;
            for (bool positive: {true, false}) {
                for (int64_t bitsOfA: {1, 2, 10, 30, 63}) {
                    BigTorusParams in_params(in_nblimbs);
                    BigTorusParams out_params(out_nblimbs);

                    BigTorus a(in_params);
                    zero(a);
                    random(a, in_nblimbs);
                    RR::SetPrecision((in_nblimbs) * BITS_PER_LIMBS);
                    RR aa = to_RR(a);

                    RR::SetOutputPrecision(100);

                    //cout << endl;
                    //cout << "aa: " << aa << endl;

                    BigTorus out(out_params);
                    //RR::SetPrecision((in_nblimbs+0) * BITS_PER_LIMBS);
                    //RR oout = (random_RR() - 0.5);
                    //to_torus(out, oout);

                    //cout << "oout: " << oout << endl;


                    //test positive first
                    ZZ coef_zz = RandomBits_ZZ(bitsOfA);
                    int64_t coef = to_long(coef_zz);
                    if (!positive) {
                        coef = -coef;
                        coef_zz = -coef_zz;
                    }

                    //cout << "coef_zz: " << coef_zz << endl;

                    mulS64(out, coef, a, limb_prec);

                    //test that out = a * coef modulo 1
                    RR::SetPrecision((in_nblimbs + 1) * BITS_PER_LIMBS);
                    RR target = centerMod1(to_RR(coef_zz) * aa);
                    RR actual = to_RR(out);

                    //cout << "actual: " << actual << endl;
                    //cout << "target: " << target << endl;

                    RR distance = centerMod1(target - actual);
                    //cout << "distance: " << distance << endl;

                    ASSERT_LE(log2DiffMod1(target, actual), -limb_prec * BITS_PER_LIMBS + bitsOfA);
                }
            }
        }
    }
}
