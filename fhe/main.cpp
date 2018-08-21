#include <iostream>
#include "commons.h"
#include "arithmetic.h"

NTL_CLIENT;

const UINT64 X_limbs = 1;
const UINT64 X_plainExp = 0;
const UINT64 X_levelExp = 20;
const BigTorusParams X_params(X_limbs, X_plainExp, X_levelExp);

const int64_t S_limbs = 1;
const int64_t S_plainExp = 1;
const int64_t S_levelExp = 20;
const BigTorusParams S_params(S_limbs, S_plainExp, S_levelExp);

const int64_t y_limbs = 1;
const int64_t y_plainExp = 1;
const int64_t y_levelExp = 20;
const BigTorusParams y_params(y_limbs, y_plainExp, y_levelExp);

const int64_t p_limbs = 1;
const int64_t p_plainExp = 1;
const int64_t p_levelExp = 20;
const BigTorusParams p_params(p_limbs, p_plainExp, p_levelExp);

const int64_t w_limbs = 1;
const int64_t w_plainExp = -1;
const int64_t w_levelExp = 20;
const BigTorusParams w_params(w_limbs, w_plainExp, w_levelExp);

const int64_t ymp_limbs = 1;
const int64_t ymp_plainExp = 1;
const int64_t ymp_levelExp = 20;
const BigTorusParams ymp_params(ymp_limbs, ymp_plainExp, ymp_levelExp);

const int64_t beta_limbs = 1;
const int64_t beta_plainExp = 6; //verify [-32,32]?
const int64_t beta_levelExp = 20;
const BigTorusParams beta_params(beta_limbs, beta_plainExp, beta_levelExp);

const int64_t XBeta_limbs = 1;
const int64_t XBeta_plainExp = 6; //verify
const int64_t XBeta_levelExp = 20;
const BigTorusParams XBeta_params(XBeta_limbs, XBeta_plainExp, XBeta_levelExp);

const int64_t mgrad_limbs = 1;
const int64_t mgrad_plainExp = 6; //verify
const int64_t mgrad_levelExp = 20;
const BigTorusParams mgrad_params(mgrad_limbs, mgrad_plainExp, mgrad_levelExp);

int main() {
    int k = 4;
    int m = 10643;
    int n = 245;
    int ITERS = 7;    //num of logreg iters
    //double step = 4.; //learning rate (close to 4)

    BigTorusMatrix X(n, k, X_params);
    BigTorusMatrix S(n, m, S_params);
    BigTorusVector y(n, y_params);
    fill_matrix_S(S);
    fill_matrix_Xy(X, y);

    cout << "y:" << y << endl;

    //
    BigTorusVector p(n, p_params);
    BigTorusVector w(n, w_params);
    BigTorusVector beta(k, beta_params);

    //logistic regression
    zero(beta);
    for (int iter = 0; iter <= ITERS; iter++) {
        cout << "--- iter " << iter << "-----" << endl;
        cout << "beta: " << beta << endl;
        BigTorusVector XBeta(n, XBeta_params);
        fixp_Ab_prod(XBeta, X, beta);
        cout << "XBeta: " << XBeta << endl;
        fixp_sigmoid_vec(p, w, XBeta);
        cout << "p: " << p << endl;
        cout << "w: " << w << endl;
        BigTorusVector ymp(n, ymp_params);
        fixp_sub(ymp, y, p);
        cout << "ymp: " << ymp << endl;
        BigTorusVector mgrad(k, mgrad_params);
        fixp_tAb_prod(mgrad, X, ymp);
        cout << "mgrad: " << mgrad << endl;
        fixp_public_scale(mgrad, 4);
        cout << "4mgrad: " << mgrad << endl;
        fixp_add(beta, beta, mgrad);
        cout << "grad: " << iter << " " << fixp_debug_norm(mgrad) << endl;
    }
    // p and W
    cout << "p: " << p << endl;
    cout << "w: " << w << endl;
    /*
    // Zstar
    vec_float zStar;
    zStar = (y - p) * S;
    cout << "zStar:" << zStar << endl;
    // G
    mat_float G = transpose(X) * diagProd(w, X);
    mat_float A = transpose(X) * diagProd(w, S);
    cout << "G: " << G << endl;
    cout << "A: " << A << endl;
    //cholesky
    for (int i = 0; i < k; i++) {
        assert(sqrt(G[i][i]) > 0); //the matrix must be positive definite!!
        float alpha = 1. / sqrt(G[i][i]);
        symetricScale(G, i, alpha);
        rowScale(A, i, alpha);
        for (int j = i + 1; j < k; j++) {
            float gamma = -G[j][i];
            symetricTransvect(G, j, i, gamma);
            rowTransvect(A, j, i, gamma);
        }
    }
    //denominator
    vec_float sStar2 = scaledColSqNorms(S, w) - colSqNorms(A);
    cout << "sStar2: " << sStar2 << endl;
    vec_float ri;
    ri.SetLength(m);
    for (int j = 0; j < m; j++)
        ri[j] = 2 * log(abs(zStar[j])) - log(abs(sStar2[j]));
    cout << "ri: " << ri << endl;
    vec_float pval;
    pval.SetLength(m);
    for (int j = 0; j < m; j++)
        pval[j] = pvalexp(ri[j]);
    cout << "pval: " << pval << endl;
    return 0;
     */
}
