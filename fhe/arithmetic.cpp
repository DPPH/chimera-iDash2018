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
#include "BigFixP.h"
#include "commons.h"

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


bool is_binary(const RR &value) {
    return value == 0 || value == 1;
}

void fill_matrix_S(BigFixPMatrix &S) {
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

void fill_matrix_Xy(BigFixPMatrix &X, BigFixPVector &y) {
    const long n = X.rows;
    const long k = X.cols;
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
    for (int i = 0; i < n; i++) {
        std::getline(ifs, line);
        for (int j = 0; j < int(line.size()); j++) {
            if (line[j] == ',') line[j] = ' ';
        }
        istringstream iss(line);
        iss >> buf; //ignore label
        iss >> rrbuf; //read y
        assert(is_binary(rrbuf)); //y binary
        to_fixP(y(i), rrbuf);

        B[0][i] = 1.; //intercept first
        nbels[0]++;
        for (int j = 1; j < k; j++) {
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
    for (int j = 0; j < k; j++) {
        sums[j] /= nbels[j];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            if (B[j][i] < 1e-75) {
                B[j][i] = sums[j];
            }
        }
    }
    //orthogonalize B
    for (int i = 0; i < k; i++) {
        //remove component on previous vectors
        for (int j = 0; j < i; j++) {
            B[i] -= (B[i] * B[j]) * B[j];
        }
        //normalize B[i]
        RR alpha = inv(sqrt(B[i] * B[i]));
        B[i] *= alpha;
    }
    //cout << B << endl;
    //cout << B*transpose(B) << endl;
    //put B in X
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            to_fixP(X(i, j), B[j][i]);
        }
    }
}

void sigmoid_vec(BigFixPVector &p, BigFixPVector &w, BigFixPVector &x) {
    assert_dramatically(p.length == x.length && p.length == w.length, "wrong dimensions");
    //this function is bootstrapped, so we will just use the RR conversion
    for (int64_t i = 0; i < p.length; i++) {
        RR xx = to_RR(x(i));
        RR sigmo = inv(1 + exp(xx));
        to_fixP(p(i), sigmo);
        to_fixP(w(i), sigmo * (1 - sigmo));
    }
}

void public_scale(BigFixPVector &res, int alpha) {
    const uint64_t length = res.length;
    const uint64_t nblimbs = res.params->torus_params.torus_limbs;
    for (uint64_t i = 0; i < length; i++) {
        bigTorusRawScale(res.limbs_raw + i * nblimbs, alpha, nblimbs);
    }
}

NTL::RR debug_norm(const BigFixPVector &v) {
    //this function is bootstrapped, so we will just use the RR conversion
    RR reps;
    reps = 0;
    for (int64_t i = 0; i < v.length; i++) {
        RR xx = to_RR(v(i));
        reps += xx * xx;
    }
    return sqrt(reps);
}

std::ostream &operator<<(std::ostream &out, const BigFixPVector &v) {
    //this function is bootstrapped, so we will just use the RR conversion
    out << "[";
    for (int64_t i = 0; i < v.length; i++) {
        out << to_double(to_RR(v(i))) << " ";
    }
    out << "]";
    return out;
}

void tAb_prod(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b) {
    const uint64_t n = A.rows;
    const uint64_t m = A.cols;
    assert_dramatically(b.length == n, "wrong dimension");
    assert_dramatically(res.length == m, "wrong dimension");
    //this function is bootstrapped, so we will just use the RR conversion
    mat_RR tAA;
    tAA.SetDims(m, n);
    vec_RR bb;
    bb.SetLength(n);
    for (int64_t i = 0; i < m; i++) {
        for (int64_t j = 0; j < n; j++) {
            tAA[i][j] = to_RR(A(j, i));
        }
    }
    for (int64_t j = 0; j < n; j++) {
        bb[j] = to_RR(b(j));
    }
    vec_RR reps = tAA * bb;
    for (int64_t i = 0; i < m; i++) {
        to_fixP(res(i), reps[i]);
    }
}
