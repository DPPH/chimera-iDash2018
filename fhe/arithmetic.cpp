#include <NTL/RR.h>
#include <NTL/ZZ_limbs.h>
#include <assert.h>
#include <cassert>
#include <gmp.h>
#include <curses.h>
#include <ncurses.h>
#include <fstream>
#include <sstream>
#include "arithmetic.h"
#include "commons.h"
#include "BigTorus.h"


NTL_CLIENT;

void to_fixP(BigTorusRef reps, const NTL::RR &a) {
    const int64_t n = reps.params.torus_limbs;
    RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    BigTorusRef ta(reps.limbs_end, reps.params);
    to_torus(ta, a * pow(2, -(reps.params.plaintext_expo + reps.params.level_expo)));
}

void to_torus(BigTorusRef reps, const NTL::RR &a) {
    assert(abs(a) <= 0.5);
    const int64_t n = reps.params.torus_limbs;
    RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    ZZ az = RoundToZZ((a + 2) * pow(2, n * BITS_PER_LIMBS));
    //cout << "az:" << az << endl;
    const UINT64 *limbs = ZZ_limbs_get(az);
    mpn_copyi(reps.limbs_end - n, limbs, n);
}

NTL::RR fixp_to_RR(const BigTorusRef &a) {
    BigTorusRef ta(a.limbs_end, a.params);
    RR reps = to_RR(ta);
    reps *= NTL::pow(to_RR(2), to_RR(long(a.params.level_expo + a.params.plaintext_expo)));
    return reps;
}

NTL::RR to_RR(const BigTorusRef &a) {
    const long n = a.params.torus_limbs;
    NTL::RR pow2m32 = NTL::pow(to_RR(2), to_RR(-32));
    static const UINT64 mask32 = 0xFFFFFFFFul;
    NTL::RR::SetPrecision(n * BITS_PER_LIMBS + 2);
    NTL::RR reps;
    reps = 0;
    for (long i = 0; i < n; i++) {
        UINT64 limb = a.limbs_end[-n + i];
        reps += (limb & mask32);
        reps *= pow2m32;
        reps += (limb >> 32u);
        reps *= pow2m32;
    }
    if (reps >= 0.5) reps -= 1;
    return reps;
}


bool is_binary(const RR &value) {
    return value == 0 || value == 1;
}

void fill_matrix_S(BigTorusMatrix &S) {
    const long m = S.cols;
    const long n = S.rows;
    basic_ifstream<char> ifs("data/snpMat.txt");
    assert_dramatically(bool(ifs), "data/snpMat.txt not found!");
    RR rrbuf;
    basic_string<char> line;
    std::getline(ifs, line); //ignore first header line
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            ifs >> rrbuf;
            assert(is_binary(rrbuf)); //S binary
            to_fixP(S(i, j), rrbuf);
        }
    assert(ifs); //verify that all values have been read
    ifs.close();
}

void fill_matrix_Xy(BigTorusMatrix &X, BigTorusVector &y) {
    const UINT64 n = X.rows;
    const UINT64 k = X.cols;
    assert_dramatically(y.length == n, "dimension of X, y are wrong");
    basic_ifstream<char> ifs("data/covariates.csv");
    assert(ifs);
    mat_RR B;
    vec_RR sums;
    vec_RR nbels;
    B.SetDims(k, n);
    sums.SetLength(k);
    clear(sums);
    nbels.SetLength(k);
    clear(nbels);

    string line;
    string buf;
    RR rrbuf;
    std::getline(ifs, line); //ignore first header line
    for (UINT64 i = 0; i < n; i++) {
        std::getline(ifs, line);
        for (int j = 0; j < int(line.size()); j++) {
            if (line[j] == ',') line[j] = ' ';
        }
        istringstream iss(line);
        iss >> buf; //ignore label
        iss >> rrbuf; //read y
        assert(is_binary(rrbuf)); //y binary
        to_fixP(y.getAT(i), rrbuf);

        B[0][i] = 1.; //intercept first
        nbels[0]++;
        for (UINT64 j = 1; j < k; j++) {
            iss >> buf;
            if (buf == "NA") {
                B[j][i] = -1e80;
            } else {
                B[j][i] = stod(buf);
                sums[j] += B[j][i];
                nbels[j]++;
            }
        }
        assert(iss);
    }
    assert(ifs);
    ifs.close();
    //replace all NANs by average
    //cout << B << endl;
    //cout << nbels << endl;
    for (UINT64 j = 0; j < k; j++) {
        sums[j] /= nbels[j];
    }
    for (UINT64 i = 0; i < n; i++) {
        for (UINT64 j = 0; j < k; j++) {
            if (B[j][i] < 1e-75) {
                B[j][i] = sums[j];
            }
        }
    }
    //orthogonalize B
    for (UINT64 i = 0; i < k; i++) {
        //remove component on previous vectors
        for (UINT64 j = 0; j < i; j++) {
            B[i] -= (B[i] * B[j]) * B[j];
        }
        //normalize B[i]
        RR alpha = inv(sqrt(B[i] * B[i]));
        B[i] *= alpha;
    }
    //cout << B << endl;
    //cout << B*transpose(B) << endl;
    //put B in X
    for (UINT64 i = 0; i < n; i++) {
        for (UINT64 j = 0; j < k; j++) {
            to_fixP(X(i, j), B[j][i]);
        }
    }
}

void fixp_sigmoid_vec(BigTorusVector &p, BigTorusVector &w, BigTorusVector &x) {
    assert_dramatically(p.length == x.length && p.length == w.length, "wrong dimensions");
    //this function is bootstrapped, so we will just use the RR conversion
    for (UINT64 i = 0; i < p.length; i++) {
        RR xx = fixp_to_RR(x.getAT(i));
        RR sigmo = inv(1 + exp(-xx));
        to_fixP(p.getAT(i), sigmo);
        to_fixP(w.getAT(i), sigmo * (1 - sigmo));
    }
}

void fixp_public_scale(BigTorusVector &res, int alpha) {
    const UINT64 length = res.length;
    const UINT64 nblimbs = res.btp.torus_limbs;
    for (UINT64 i = 0; i < length; i++) {
        bigTorusRawScaleU64(res.limbs_end - nblimbs + i * nblimbs, alpha, nblimbs);
    }
}

NTL::RR fixp_debug_norm(const BigTorusVector &v) {
    //this function is bootstrapped, so we will just use the RR conversion
    RR reps;
    reps = 0;
    for (UINT64 i = 0; i < v.length; i++) {
        RR xx = fixp_to_RR(v.getAT(i));
        reps += xx * xx;
    }
    return sqrt(reps);
}

std::ostream &operator<<(std::ostream &out, const BigTorusVector &v) {
    //this function is bootstrapped, so we will just use the RR conversion
    out << "[";
    for (UINT64 i = 0; i < v.length; i++) {
        out << to_double(fixp_to_RR(v.getAT(i))) << " ";
        if (i == 5 && v.length > 10) {
            out << "... ";
            i = v.length - 6;
        }
    }
    out << "]";
    return out;
}

void fixp_tAb_prod(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b) {
    fixp_tAb_prod_fake(res, A, b);
}

void fixp_tAb_prod_fake(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b) {
    const UINT64 n = A.rows;
    const UINT64 m = A.cols;
    assert_dramatically(b.length == n, "wrong dimension");
    assert_dramatically(res.length == m, "wrong dimension");
    //this function is bootstrapped, so we will just use the RR conversion
    mat_RR tAA;
    tAA.SetDims(m, n);
    vec_RR bb;
    bb.SetLength(n);
    for (UINT64 i = 0; i < m; i++) {
        for (UINT64 j = 0; j < n; j++) {
            tAA[i][j] = fixp_to_RR(A(j, i));
        }
    }
    for (UINT64 j = 0; j < n; j++) {
        bb[j] = fixp_to_RR(b.getAT(j));
    }
    vec_RR reps = tAA * bb;
    for (UINT64 i = 0; i < m; i++) {
        to_fixP(res.getAT(i), reps[i]);
    }
}

void fixp_Ab_prod(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b) {
    fixp_Ab_prod_fake(res, A, b);
}

void fixp_Ab_prod_fake(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b) {
    const UINT64 n = A.rows;
    const UINT64 m = A.cols;
    assert_dramatically(b.length == m, "wrong dimension");
    assert_dramatically(res.length == n, "wrong dimension");
    //this function is bootstrapped, so we will just use the RR conversion
    mat_RR AA;
    AA.SetDims(n, m);
    vec_RR bb;
    bb.SetLength(m);
    for (UINT64 i = 0; i < n; i++) {
        for (UINT64 j = 0; j < m; j++) {
            AA[i][j] = fixp_to_RR(A(i, j));
        }
    }
    for (UINT64 j = 0; j < m; j++) {
        bb[j] = fixp_to_RR(b.getAT(j));
    }
    vec_RR reps = AA * bb;
    for (UINT64 i = 0; i < n; i++) {
        to_fixP(res.getAT(i), reps[i]);
    }
}

uint8_t random_bit() {
    return uint8_t(random_uint64_t() % 2);
}

UINT64 random_uint64_t() {
    UINT64 buf;
    mpn_random(&buf, 1);
    return buf;
}
