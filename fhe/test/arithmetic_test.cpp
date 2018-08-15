#include "gtest/gtest.h"
#include <NTL/LLL.h>
#include "../commons.h"
#include "../arithmetic.h"
#include "../BigReal.h"
#include "../BigComplex.h"
#include "../TRLwe.h"

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

TEST(BIGTORUS_ARITHMETIC, submul64) {
    for (int64_t limb_prec:  {1, 3, 5, 7}) {
        for (int64_t out_nblimbs: {limb_prec, limb_prec + 1, limb_prec + 3}) {
            for (int64_t in_nblimbs: {limb_prec + 1, limb_prec + 2, limb_prec + 3}) {
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
                    ZZ coef_zz = RandomBits_ZZ(32);
                    int64_t coef = to_long(coef_zz);
                    if (!positive) {
                        coef = -coef;
                        coef_zz = -coef_zz;
                    }

                    //cout << "coef_zz: " << coef_zz << endl;

                    subMulS64(out, coef, a, limb_prec);

                    //test that out = oldout - a * coef modulo 1
                    RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
                    RR target = centerMod1(oout - to_RR(coef_zz) * aa);
                    RR actual = to_RR(out);

                    //cout << "actual: " << actual << endl;
                    //cout << "target: " << target << endl;

                    RR distance = centerMod1(target - actual);
                    //cout << "distance: " << distance << endl;

                    ASSERT_LE(log2DiffMod1(target, actual), -limb_prec * BITS_PER_LIMBS + 1 + 32);
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

TEST(BIGTORUS_ARITHMETIC, BigReal_to_BigTorus) {
    for (const int64_t bits_a : {1, 17, 31}) {
        for (const int64_t nblimbs : {1, 2, 5}) {
            for (const int64_t rlimbs: {nblimbs, nblimbs + 1}) {
                BigTorusParams params(nblimbs);

                BigReal source(rlimbs);
                BigTorus dest(params);

                RR::SetPrecision(rlimbs * BITS_PER_LIMBS);
                RR rsource = random_RR();
                to_BigReal(source, rsource);
                // expected: source * 2^a centermod 1
                RR rbi = to_RR(source) * pow(to_RR(2), to_RR(long(bits_a)));
                rbi -= to_RR(RoundToZZ(rbi));
                // actual:
                to_BigTorus(dest, source, bits_a, nblimbs);
                RR outi = to_RR(dest);
                //RR::SetOutputPrecision(nblimbs * BITS_PER_LIMBS / log2(10.));
                //cout << rbi << endl;
                //cout << outi << endl;
                ASSERT_LE(log2Diff(rbi, outi), -nblimbs * BITS_PER_LIMBS);
            }
        }
    }
}

TEST(BIGCOMPLEX_ARITHMETIC, mulTo) {
    for (const int64_t nblimbs : {1, 2, 5}) {
        RR::SetPrecision(nblimbs * BITS_PER_LIMBS + 10);
        RR care = random_RR();
        RR caim = random_RR();
        RR cbre = random_RR();
        RR cbim = random_RR();
        BigComplex ca(nblimbs);
        BigComplex cb(nblimbs);
        to_BigReal(ca.real, care);
        to_BigReal(cb.real, cbre);
        to_BigReal(ca.imag, caim);
        to_BigReal(cb.imag, cbim);
        mulTo(cb, ca);
        RR targetre = cbre * care - cbim * caim;
        RR targetim = cbre * caim + cbim * care;
        RR actualre = to_RR(cb.real);
        RR actualim = to_RR(cb.imag);
        ASSERT_LE(log2Diff(targetre, actualre), -nblimbs * BITS_PER_LIMBS + 1);
        ASSERT_LE(log2Diff(targetim, actualim), -nblimbs * BITS_PER_LIMBS + 2);
    }
}

TEST(BIGTORUS_ARITHMETIC, BigTorus_To_BigReal) {
    for (const int64_t nblimbs : {1, 2, 5}) {
        RR::SetPrecision(nblimbs * BITS_PER_LIMBS + 10);

        BigTorusParams params(nblimbs);
        BigTorus b(params);
        random(b);
        // expected:
        RR bi = to_RR(b);
        // actual:
        BigReal rb(nblimbs);
        to_BigReal(rb, b);
        RR rbi = to_RR(rb);
        ASSERT_LE(log2Diff(bi, rbi), -nblimbs * BITS_PER_LIMBS);
    }
}

TEST(BIGTORUS_ARITHMETIC, int2_To_BigReal) {
    for (const int64_t bits_a : {1, 17, 31}) {
        for (const int64_t nblimbs : {1, 2, 5}) {
            BigReal dest(nblimbs);

            int64_t bit_mask = (1ul << bits_a) - 1;
            int64_t bit_offset = (1ul << (bits_a - 1)) - 1;
            int64_t source = (random() & bit_mask) - bit_offset;

            RR::SetPrecision(nblimbs * BITS_PER_LIMBS + 10);
            RR ai = to_RR(long(source)) / pow(to_RR(2), to_RR(long(bits_a)));
            to_BigReal(dest, source, bits_a);
            RR rai = to_RR(dest);
            ASSERT_LE(log2Diff(ai, rai), -nblimbs * BITS_PER_LIMBS);
        }
    }
}

TEST(BIGTORUS_ARITHMETIC, bitdecomp32) {
    for (const int64_t nblimbs : {1, 2, 5}) {
        BigTorusParams params(nblimbs);
        int ell = 2 * nblimbs;
        BigTorus b(params);
        BigTorus bof(params);
        BigTorus hi(params);
        BigTorus recomp_b(params);
        BigTorus tmp(params);
        random(b);

        zero(recomp_b);
        bitdecomp_signed_offset32_apply(bof, b);
        for (int i = 0; i < ell; i++) {
            int64_t coeff = bitdecomp_signed_coef32(bof, i + 1, nblimbs);
            zero(hi);
            hi.limbs_end[-(i + 2) / 2] = (i % 2 == 0) ? 0x100000000ul : 0x1ul;
            mulS64(tmp, coeff, hi, nblimbs);
            add(recomp_b, recomp_b, tmp);
        }
        ASSERT_LE(log2Diff(b, recomp_b), -nblimbs * BITS_PER_LIMBS);
    }
}
