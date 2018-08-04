#include "gtest/gtest.h"
#include <NTL/LLL.h>
#include "../BigFixP.h"
#include "../commons.h"
#include "../arithmetic.h"

NTL_CLIENT;

double log2Diff(const RR &a, const RR &b) {
    if (a == b) return -INFINITY;
    return to_double(log(abs(a - b))) / log(2.);
}

TEST(BIGTORUS_ARITHMETIC, convRR) {
    for (int torus_limbs: {5, 10, 15}) {
        BigTorusParams params(torus_limbs);

        //torus load and store to RR
        for (int trial = 0; trial < 100; trial++) {
            RR test1 = random_RR() - 0.5;
            BigTorus testTorus(&params);
            to_torus(testTorus, test1);
            RR test2 = to_RR(testTorus);

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
                BigFixPParams params(torus_limbs, plaintext_expo, level_expo);
                int64_t torus_prec = torus_limbs * BITS_PER_LIMBS - plaintext_expo - level_expo;

                //fixP load and store to RR
                for (int trial = 0; trial < 100; trial++) {
                    RR test1 = 1000 * (random_RR() - 0.5);
                    BigFixP testFixp(&params);
                    to_fixP(testFixp, test1);
                    RR test2 = to_RR(testFixp);

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
    BigFixPParams params(10, 10, 128);
    BigFixPParams paramsb(11, 10, 128);
    BigFixPParams paramsc(10, 10, 64);

    //fixP addition
    for (int trial = 0; trial < 100; trial++) {
        RR testa = 1000 * (random_RR() - 0.5);
        RR testb = 1000 * (random_RR() - 0.5);
        BigFixP testFixpA(&params);
        BigFixP testFixpB(&paramsb);
        BigFixP testFixpC(&paramsc);
        to_fixP(testFixpA, testa);
        to_fixP(testFixpB, testb);
        add(testFixpC, testFixpA, testFixpB, 640);
        RR testc = to_RR(testFixpC);

        RR::SetOutputPrecision(long(640 * log(2) / log(10)));
        //cout << testa << endl;
        //cout << testb << endl;
        //cout << testa+testb << endl;
        //cout << testc << endl;
        ASSERT_LE(log2Diff(testa + testb, testc), -500);
    }
}
