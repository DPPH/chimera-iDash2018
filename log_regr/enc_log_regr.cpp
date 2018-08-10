#include "io.h"
#include "common.h"
#include "lr_params.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"

#include <vector>

using namespace std;

/**
 * @brief Build a test-polynomial for X[row,j] * \sigmoid(u) j=0..k-1
 */
void build_test_poly_sigmoid_xt(TorusPolynomial<Torus>* test_poly, const Data& data, const LRParams& lr_params, int row) {
  const auto& Xi = data.X[row];
  const int& k = data.k;

  const int& sigm_levels = lr_params.sigm_levels;
  const int& sigm_min = lr_params.sigm_in_min;
  const int& sigm_max = lr_params.sigm_in_max;
  const double& scale = lr_params.scale;

  int N = test_poly->N;
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

void encrypt_sigmoid_xt(TLweSample<Torus> *test_poly_encs, const Data& data, const LRParams& lr_params, const TLweParams<Torus>* params, const TLweKey<Torus>* key) {
  int N = params->N;
  double alpha = lr_params.alpha_test_poly;

  TorusPolynomial<Torus> *test_poly = new_obj<TorusPolynomial<Torus>>(N);
  for (int i = 0; i < data.n; ++i) {
    TorusPolyFunctions<Torus>::Clear(test_poly);
    build_test_poly_sigmoid_xt(test_poly, data, lr_params, i);

    TLweFunctions<Torus>::SymEncrypt(test_poly_encs+i, test_poly, alpha, key);
  }

  del_obj(test_poly);
}

/**
 * @brief Encrypt y vector as \sum_i y_i Z^(n-i)
 */
void encrypt_y(TLweSample<Torus> *enc_y, const Data& data, const LRParams& lr_params, const TLweParams<Torus>* params, const TLweKey<Torus>* key) {
    int N = params->N;
    double alpha = lr_params.alpha_test_poly; //TODO: change to good one

    assert(data.n < N);

    TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
    Torus *coefs = msg->coefsT;
    for (int i = 0; i < data.n; ++i) {
        coefs[i] = data.y[data.n-i];
    }

    TLweFunctions<Torus>::SymEncrypt(enc_y, msg, alpha, key);

    del_obj(msg);
}

/**
 * @brief Encrypt X columns as \sum_i X_ij Z^i for j=0..k-1
 */
void encrypt_X_cols(TGswSample<Torus>* enc_X_cols, const Data& data, const LRParams& lr_params, const TGswParams<Torus>* params, const TGswKey<Torus>* key) {
    int N = params->tlwe_params->N;
    double alpha = lr_params.alpha_test_poly; //TODO: change to good one

    assert(data.n < N);

    IntPolynomial *msg = new_obj<IntPolynomial>(N);
    int *coefs = msg->coefs;

    for (int j = 0; j < data.k; ++j) {
        for (int i = 0; i < data.n; ++i) {
            coefs[i] = data.X[i][j];
        }
        TGswFunctions<Torus>::SymEncrypt(enc_X_cols+j, msg, alpha, key);
    }

    del_obj(msg);
}

int main(int argc, char const *argv[]) {
    Data data;
    fill_data(data);

    LRParams lr_params;
    lr_params.n = data.n;
    lr_params.k = data.k;
    lr_params.m = data.m;

    int N = 8192>>1;
    TorusPolynomial<Torus> *test_poly = new_TorusPolynomial<Torus>(N);
    build_test_poly_sigmoid_xt(test_poly, data, lr_params, 1);
}

