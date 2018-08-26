#include <include/gtest/gtest.h>

#include "../log_regr_fncs.h"
#include "../tlwekeyswitch.h"
#include "../tlwe-functions-extra.h"
#include "test_common.h"
#include <NTL/ZZX.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

class LogRegrFuncTest: public TestCommon {

};


// TEST_F(LogRegrFuncTest, compute_Xt_y) {
//     const TGswKey<Torus>* trgsw_key = secret_keyset->trgsw_key_l2;
//     const TLweKey<Torus>* trlwe_key = secret_keyset->trlwe_key_l2;;
//     const LweKey<Torus>* tlwe_key = secret_keyset->tlwe_key_l2;

//     const TGswParams<Torus>* trgsw_params = trgsw_key->params;
//     const TLweParams<Torus>* trlwe_params = trlwe_key->params;
//     const LweParams<Torus>* tlwe_params = tlwe_key->params;

//     const int k = 4;
//     const int n = 16;
//     lr_params.k = k;
//     lr_params.n = n;

//     /* encrypt X_cols */
//     Vec<ZZX> clear_X_cols;
//     TGswSample<Torus>* X_cols = nullptr;
//     {
//         clear_X_cols.SetLength(k);
//         X_cols = new_obj_array<TGswSample<Torus>>(k, trgsw_params);

//         cout << "X_cols:" << endl;
//         for (int i = 0; i < k; ++i) {
//             IntPolynomial* msg = new_obj<IntPolynomial>(trlwe_params->N);
//             for (int j = 0; j < n; ++j) {
//                 msg->coefs[j] = i*n+j + 1;
//                 SetCoeff(clear_X_cols[i], j, msg->coefs[j]);
//             }
//             TGswFunctions<Torus>::SymEncrypt(X_cols+i, msg, trlwe_params->alpha_min, trgsw_key);
//             del_obj(msg);

//             cout << clear_X_cols[i] << endl;
//         }
//     }

//     /* encrypt y */
//     ZZX clear_y;
//     TLweSample<Torus>* y = nullptr;
//     {
//         TorusPolynomial<Torus>* tmp = new_obj<TorusPolynomial<Torus>>(trlwe_params->N);
//         for (int j = 0; j < n; ++j)  {
//             // tmp->coefsT[j] = j + 1;
//             tmp->coefsT[j] = TorusUtils<Torus>::from_double(j / 4. / n);
//             SetCoeff(clear_y, j, tmp->coefsT[j]);
//         }

//         y = new_obj<TLweSample<Torus>>(trlwe_params);
//         TLweFunctions<Torus>::SymEncrypt(y, tmp, trlwe_params->alpha_min, trlwe_key);
//         cout << "y:" << endl << clear_y << endl;
//     }

//     /* perform computation */
//     LweSample<Torus>* Xt_y = new_obj_array<LweSample<Torus>>(k, tlwe_params);
//     compute_Xt_y(Xt_y, X_cols, y, trgsw_params, lr_params);

//     /* verify result */
//     for (int i = 0; i < k; ++i) {

//         Torus t = LweFunctions<Torus>::SymDecrypt(Xt_y+i, tlwe_key, 4 * n * 256);
//         double comp_res = TorusUtils<Torus>::to_double(t);
//         printf("%d %lf %ld\n", i, comp_res, t);

//         ZZX tmp = clear_X_cols[i] * clear_y;
//         // cout << tmp << endl;
//         cout << tmp[n-1] << endl;
//         cout << tmp[n-1] % (1L<<63) << endl;
//         // cout << TorusUtils<Torus>::to_double(tmp[n-1] % (1L<<63)) << endl;
//     }

//     del_obj_array(k, Xt_y);
//     del_obj(y);
//     del_obj_array(k, X_cols);
// }


// TEST_F(LogRegrFuncTest, ks_l1_l0) {
//     const LweKey<Torus>* tlwe_key_l0 = secret_keyset->tlwe_key_l0;
//     const LweKey<Torus>* tlwe_key_l1 = secret_keyset->tlwe_key_l1;

//     const LweParams<Torus>* tlwe_params_l0 = tlwe_key_l0->params;
//     const LweParams<Torus>* tlwe_params_l1 = tlwe_key_l1->params;

//     RandomGen::set_seed(41);

//     // LweSample<Torus>* inp_sample1 = read_tlwe_samples("ks_inp.ctxt", tlwe_params_l1, 1);

//     LweSample<Torus>* inp_sample2 = new_obj<LweSample<Torus>>(tlwe_params_l1);
//     double msg = 0.5911391 * lr_params.X_beta_scale;
//     LweFunctions<Torus>::SymEncrypt(inp_sample2, TorusUtils<Torus>::from_double(msg), pow(2.0,-16), tlwe_key_l1);

//     LweSample<Torus>* inp_sample = inp_sample2;

//     printf("inp tlwe sample: ");
//     print_tlwe_sample(inp_sample, tlwe_key_l1, 1./lr_params.X_beta_scale);
//     printf("\n");

//     LweSample<Torus>* out_sample = new_obj<LweSample<Torus>>(tlwe_params_l0);
//     LweFunctions<Torus>::KeySwitch(out_sample, cloud_keyset->ks_l1_l0, inp_sample);

//     printf("out tlwe sample: ");
//     print_tlwe_sample(out_sample, tlwe_key_l0, 1./lr_params.X_beta_scale);
//     printf("\n");

// }

TEST_F(LogRegrFuncTest, ks_l2_l1) {
    const TLweKey<Torus>* trlwe_key_l1 = secret_keyset->trlwe_key_l1;
    const LweKey<Torus>* tlwe_key_l2 = secret_keyset->tlwe_key_l2;

    const TLweParams<Torus>* trlwe_params_l1 = trlwe_key_l1->params;
    const LweParams<Torus>* tlwe_params_l2 = tlwe_key_l2->params;

    RandomGen::set_seed(41);

    // LweSample<Torus>* inp_sample1 = read_tlwe_samples("ks_inp.ctxt", tlwe_params_l1, 1);

    LweSample<Torus>* inp_sample2 = new_obj<LweSample<Torus>>(tlwe_params_l2);
    double msg = 0.1911391;
    LweFunctions<Torus>::SymEncrypt(inp_sample2, TorusUtils<Torus>::from_double(msg), pow(2.0,-20), tlwe_key_l2);

    LweSample<Torus>* inp_sample = inp_sample2;

    printf("inp tlwe sample: ");
    print_tlwe_sample(inp_sample, tlwe_key_l2, 1.);
    printf("\n");

    TLweSample<Torus>* out_sample = new_obj<TLweSample<Torus>>(trlwe_params_l1);
    TLweFunctionsExtra<Torus>::KeySwitch(out_sample, cloud_keyset->ks_l2_l1, inp_sample);

    printf("out tlwe sample: ");
    print_trlwe_sample(out_sample, trlwe_key_l1, 10, 1.);
    printf("\n");

}




// TEST_F(LogRegrFuncTest, blindrotate) {
//     const TGswKey<Torus>* trgsw_key = secret_keyset->trgsw_key_l2;
//     const TGswParams<Torus>* trgsw_params = trgsw_key->params;

//     const TLweKey<Torus>* trlwe_key = secret_keyset->trlwe_key_l2;
//     const TLweParams<Torus>* trlwe_params = trlwe_key->params;

//     const LweKey<Torus>* tlwe_key = secret_keyset->tlwe_key_l0;
//     const LweParams<Torus>* tlwe_params = tlwe_key->params;


//     LweSample<Torus>* inp = new_obj<LweSample<Torus>>(tlwe_params);

//     const int N = trlwe_params->N;
//     TorusPolynomial<Torus>* tp = new_obj<TorusPolynomial<Torus>>(N);
//     TorusPolynomial<Torus>* tp_copy = new_obj<TorusPolynomial<Torus>>(N);
//     // for (int i = 0; i < N/2; ++i) {
//     //     tp->coefsT[i] = TorusUtils<Torus>::from_double(double(N/2+i)/N/2);
//     // }
//     // for (int i = 0; i < N/2; ++i) {
//     //     tp->coefsT[N/2+i] = TorusUtils<Torus>::from_double(double(-i)/N/2);
//     // }

//     const int sigm_levels = N;
//     const int sigm_min = lr_params.sigm_in_min;
//     const int sigm_max = lr_params.sigm_in_max;

//     float dv = 1.0 / (sigm_levels-1) * (sigm_max - sigm_min);

//     /*0..1/2*/
//     // const double scale = 1./2;
//     // printf("sigmoid:\n");
//     // for (int q = 0; q < sigm_levels/2; ++q) {
//     //     float sigm_val = sigmoid((float)q * dv + sigm_min);
//     //     tp->coefsT[N/2+q] = TorusUtils<Torus>::from_double(-sigm_val * scale);
//     // }
//     // for (int q = 0; q < sigm_levels/2; ++q) {
//     //     float sigm_val = sigmoid((float)(q+sigm_levels/2) * dv + sigm_min);
//     //     tp->coefsT[q] = TorusUtils<Torus>::from_double(sigm_val * scale);
//     // }
//     // printf("\n");


//     const double scale = 1./2;
//     printf("sigmoid:\n");
//     for (int q = 0; q < sigm_levels; ++q) {
//         float sigm_val = sigmoid((float)q * dv + sigm_min);
//         tp->coefsT[q] = TorusUtils<Torus>::from_double(sigm_val * scale);
//     }
//     printf("\n");
//     TorusPolyFunctions<Torus>::Copy(tp_copy, tp);
//     TorusPolyFunctions<Torus>::MulByXai(tp, 2*N-N/2-16, tp_copy);


//     TLweSample<Torus> *acc = new_obj<TLweSample<Torus>>(trlwe_params);
//     for (double msg = -0.25; msg <= 0.25; msg+=0.05)
//     {
//         LweFunctions<Torus>::SymEncrypt(inp, TorusUtils<Torus>::from_double(msg), tlwe_params->alpha_min, tlwe_key);
//         TLweFunctions<Torus>::SymEncrypt(acc, tp, trlwe_params->alpha_min, trlwe_key);

//         blindrotate(acc, trlwe_params, inp, tlwe_params, cloud_keyset->bk_fft, trgsw_params, lr_params);

//         printf("%8.6lf -> ", msg);
//         print_trlwe_sample(acc, trlwe_key, 1);
//         printf("\n");
//     }

//     del_obj(acc);
//     del_obj(tp);
//     del_obj(inp);
// }
