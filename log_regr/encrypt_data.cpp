#include "io.h"
#include "common.h"
#include "lr_params.h"
#include "keyset.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "tfhe_io.h"

#include <iostream>
#include <fstream>

using namespace std;

/**
 * @brief Build a test-polynomial for X[row,j] * \sigmoid(u) j=0..k-1
 */
void build_test_poly_sigmoid_xt(TorusPolynomial<Torus>* test_poly, const vec_float& Xi, const LRParams& lr_params) {
  const int k = Xi.length();

  const int sigm_levels = lr_params.sigm_levels;
  const int sigm_min = lr_params.sigm_in_min;
  const int sigm_max = lr_params.sigm_in_max;
  const double scale = lr_params.scale;

  const int N = test_poly->N;
  Torus* test_poly_coefs = test_poly->coefsT;

  int d = (int)((float)N  / k / sigm_levels);
  float dv = 1.0 / (sigm_levels-1) * (sigm_max - sigm_min);
  for (int q = 0; q < sigm_levels; ++q) {
    float sigm_val = sigmoid((float)q * dv + sigm_min);

    for (int p = 0; p < k; ++p) {
      int st = (q*k+p)*d - d/2; // -d/2 to wrap arround X^0
      for (int i = st; i < st+d; ++i) {
        if (i >= 0) {
          test_poly_coefs[i] = TorusUtils<Torus>::from_double(Xi[p] * sigm_val * scale);
        } else {
          test_poly_coefs[N+i] = TorusUtils<Torus>::from_double(- Xi[p] * sigm_val * scale);
        }
      }
    }
  }
  // for (int i = 0; i < N; ++i) {
  //   printf("%d %lf\n", i, TorusUtils<Torus>::to_double(test_poly_coefs[i]));
  // }
}

/**
 * @brief Build, encrypt and write test-polynomial Xt . sigmoid
 */
void encrypt_write_sigmoid_xt(Ostream& out_stream, const Data& data, const LRParams& lr_params, const TLweKey<Torus> *key) {
    const TLweParams<Torus> *params = key->params;
    const int N = params->N;
    const double alpha = params->alpha_min;

    #pragma omp parallel for ordered schedule(static,1)
    for (int i = 0; i < data.n; ++i) {
        TorusPolynomial<Torus> *test_poly = new_obj<TorusPolynomial<Torus>>(N);
        TLweSample<Torus> *test_poly_enc = new_obj<TLweSample<Torus>>(params);

        TorusPolyFunctions<Torus>::Clear(test_poly);
        build_test_poly_sigmoid_xt(test_poly, data.X[i], lr_params);
        TLweFunctions<Torus>::SymEncrypt(test_poly_enc, test_poly, alpha, key);

        #pragma omp ordered
        IOFunctions<Torus>::write_tLweSample(out_stream, test_poly_enc, params);

        del_obj(test_poly_enc);

        #ifndef NDEBUG
            if (i==0) {
                printf("sigmoid_xt[%d]:\n", i);
                for (int j = 0; j < N; ++j)
                    printf("%lf ", TorusUtils<Torus>::to_double(test_poly->coefsT[j]));
                printf("\n");
            }
        #endif

        del_obj(test_poly);
    }

}

/**
 * @brief Encrypt and write y vector as \sum_i y_i Z^(n-i)
 */
void encrypt_write_y(Ostream& out_stream, const vec_float& y, const LRParams& lr_params, const TLweKey<Torus> *key) {
    const TLweParams<Torus> *params = key->params;
    const int N = params->N;
    const double alpha = params->alpha_min;

    const double scale = lr_params.scale;

    const int n = y.length();
    assert(n < N);

    TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
    TorusPolyFunctions<Torus>::Clear(msg);
    for (int i = 0; i < n; ++i) {
        msg->coefsT[i] = TorusUtils<Torus>::from_double(y[n-i-1] * scale);
    }

    #ifndef NDEBUG
        printf("y:\n");
        for (int j = 0; j < n; ++j)
            printf("%lf ", TorusUtils<Torus>::to_double(msg->coefsT[j]));
        printf("\n");
    #endif

    TLweSample<Torus> *enc_y = new_obj<TLweSample<Torus>>(params);
    TLweFunctions<Torus>::SymEncrypt(enc_y, msg, alpha, key);
    IOFunctions<Torus>::write_tLweSample(out_stream, enc_y, params);

    del_obj(enc_y);
    del_obj(msg);
}

/**
 * @brief Encrypt X columns as \sum_i X_ij Z^i for j=0..k-1
 */
void encrypt_write_X_cols(Ostream& out_stream, const Data& data, const LRParams& lr_params, const TGswKey<Torus> *key) {
    const TGswParams<Torus> *params = key->params;
    const int N = params->tlwe_params->N;
    const double alpha = params->tlwe_params->alpha_min;

    assert(data.n < N);

    const int prec_mult = 1<<lr_params.X_precbit;

    #pragma omp parallel for ordered schedule(static,1)
    for (int j = 0; j < data.k; ++j) {
        IntPolynomial *msg = new_obj<IntPolynomial>(N);
        IntPolyFunctions::Clear(msg);
        for (int i = 0; i < data.n; ++i) {
            int c = round((data.X[i][j] + 1.0) * prec_mult) - prec_mult/2;
            msg->coefs[i] = c;
        }

        TGswSample<Torus>* enc_X_col = new_obj<TGswSample<Torus>>(params);
        TGswFunctions<Torus>::SymEncrypt(enc_X_col, msg, alpha, key);

        #pragma omp ordered
        IOFunctions<Torus>::write_tGswSample(out_stream, enc_X_col, params);

        del_obj(enc_X_col);

        #ifndef NDEBUG
            printf("Xt[%d] (rescaled to int):\n", j);
            for (int i = 0; i < data.n; ++i)
                printf("%d ", msg->coefs[i]);
            printf("\n");
        #endif

        del_obj(msg);
    }
}

int main(int argc, char const *argv[]) {
    Data data;
    fill_data(data);

    LRParams lr_params;
    lr_params.n = data.n;
    lr_params.k = data.k;
    lr_params.m = data.m;

    RandomGen::set_seed(42);

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);

    ofstream out(lr_params.filename_data, ofstream::binary);
    StdOstream out_stream = to_Ostream(out);

    out_stream.fwrite(&data.n, sizeof(int32_t));
    out_stream.fwrite(&data.m, sizeof(int32_t));
    out_stream.fwrite(&data.k, sizeof(int32_t));

    encrypt_write_sigmoid_xt(out_stream, data, lr_params, keyset->trlwe_key_l2);
    encrypt_write_y(out_stream, data.y, lr_params, keyset->trlwe_key_l2);
    encrypt_write_X_cols(out_stream, data, lr_params, keyset->trgsw_key_l1);
    encrypt_write_X_cols(out_stream, data, lr_params, keyset->trgsw_key_l2);

    return 0;
}

