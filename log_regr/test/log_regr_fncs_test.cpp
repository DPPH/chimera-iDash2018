#include <include/gtest/gtest.h>

#include "../log_regr_fncs.h"
#include "test_common.h"
#include <NTL/ZZX.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

class LogRegrFuncTest: public TestCommon {

};


TEST_F(LogRegrFuncTest, compute_Xt_y) {
    const TGswKey<Torus>* trgsw_key = secret_keyset->trgsw_key_l2;
    const TLweKey<Torus>* trlwe_key = secret_keyset->trlwe_key_l2;;
    const LweKey<Torus>* tlwe_key = secret_keyset->tlwe_key_l2;

    const TGswParams<Torus>* trgsw_params = trgsw_key->params;
    const TLweParams<Torus>* trlwe_params = trlwe_key->params;
    const LweParams<Torus>* tlwe_params = tlwe_key->params;

    const int k = 4;
    const int n = 16;
    lr_params.k = k;
    lr_params.n = n;

    /* encrypt X_cols */
    Vec<ZZX> clear_X_cols;
    TGswSample<Torus>* X_cols = nullptr;
    {
        clear_X_cols.SetLength(k);
        X_cols = new_obj_array<TGswSample<Torus>>(k, trgsw_params);

        cout << "X_cols:" << endl;
        for (int i = 0; i < k; ++i) {
            IntPolynomial* msg = new_obj<IntPolynomial>(trlwe_params->N);
            for (int j = 0; j < n; ++j) {
                msg->coefs[j] = i*n+j + 1;
                SetCoeff(clear_X_cols[i], j, msg->coefs[j]);
            }
            TGswFunctions<Torus>::SymEncrypt(X_cols+i, msg, trlwe_params->alpha_min, trgsw_key);
            del_obj(msg);

            cout << clear_X_cols[i] << endl;
        }
    }

    /* encrypt y */
    ZZX clear_y;
    TLweSample<Torus>* y = nullptr;
    {
        TorusPolynomial<Torus>* tmp = new_obj<TorusPolynomial<Torus>>(trlwe_params->N);
        for (int j = 0; j < n; ++j)  {
            // tmp->coefsT[j] = j + 1;
            tmp->coefsT[j] = TorusUtils<Torus>::from_double(j / 4. / n);
            SetCoeff(clear_y, j, tmp->coefsT[j]);
        }

        y = new_obj<TLweSample<Torus>>(trlwe_params);
        TLweFunctions<Torus>::SymEncrypt(y, tmp, trlwe_params->alpha_min, trlwe_key);
        cout << "y:" << endl << clear_y << endl;
    }

    /* perform computation */
    LweSample<Torus>* Xt_y = new_obj_array<LweSample<Torus>>(k, tlwe_params);
    compute_Xt_y(Xt_y, X_cols, y, trgsw_params, lr_params);

    /* verify result */
    for (int i = 0; i < k; ++i) {

        Torus t = LweFunctions<Torus>::SymDecrypt(Xt_y+i, tlwe_key, 4 * n * 256);
        double comp_res = TorusUtils<Torus>::to_double(t);
        printf("%d %lf %ld\n", i, comp_res, t);

        ZZX tmp = clear_X_cols[i] * clear_y;
        // cout << tmp << endl;
        cout << tmp[n-1] << endl;
        cout << tmp[n-1] % (1L<<63) << endl;
        // cout << TorusUtils<Torus>::to_double(tmp[n-1] % (1L<<63)) << endl;
    }

    del_obj_array(k, Xt_y);
    del_obj(y);
    del_obj_array(k, X_cols);
}

// /**
//  * @brief Compute Xt.y from TRGSW samples of X columns and TRLWE sample of y (coefficient packed)
//  *
//  * @param Xt_y k output TLWE samples under extracted TLWE key from TRGSW sample
//  * @param X_cols k TRGSW samples encoding \sum_{i=0..n-1} Z^i . X_ij for j=0..k
//  * @param y TRLWE sample encodin \sum_{i=0..n-1} Z^(n-i) . y_i
//  * @param trgsw_params parameters of TRGSW sample
//  * @param n number of rows in matrix X
//  * @param k number of columns in matrix X
//  */
// void compute_Xt_y(
//     LweSample<Torus>* Xt_y,
//     const TGswSample<Torus>* X_cols,
//     const TLweSample<Torus>* y,
//     const TGswParams<Torus>* trgsw_params,
//     const LRParams& lr_params
//     )//
