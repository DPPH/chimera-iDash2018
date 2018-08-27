#include <cstdint>
#include <NTL/mat_RR.h>
#include <memory>
#include "TLwe.h"
#include "TRLwe.h"
#include "mainalgo.h"

NTL_CLIENT;

NTL::mat_RR debug_S;
NTL::mat_RR debug_X;
NTL::vec_RR debug_W;
shared_ptr<TLweKey> debug_key;


int main() {
    int64_t k = 3;
    //int64_t n = 13;
    int64_t n = 250;
    int64_t l = 5;
    int64_t N = 128;
    //int64_t N = 4096;

    int64_t tau_X = -5;
    int64_t tau_W = -1;
    int64_t tau_S = 1;

    int64_t L_X = 36;
    int64_t L_W = 54;
    int64_t L_S = 54;
    int64_t rho = 18;

    mat_RR plaintext_X;
    plaintext_X.SetDims(n, k + 1);
    mat_RR plaintext_S;
    plaintext_S.SetDims(n, (l * N) / 2);
    vec_RR plaintext_W;
    plaintext_W.SetLength(n);

    for (int i = 0; i < n; i++) {
        plaintext_W[i] = random_RR() * power2_RR(tau_W);
        for (int j = 0; j < k + 1; j++) {
            plaintext_X[i][j] = random_RR() * power2_RR(tau_X);
        }
        for (int j = 0; j < l * N / 2; j++) {
            plaintext_S[i][j] = random_RR() * power2_RR(tau_S);
        }
    }
    debug_S = plaintext_S;
    debug_W = plaintext_W;
    debug_X = plaintext_X;

    BigTorusParams bt_params_key(0, 0, 0);
    TRLweParams trlweParams_key(N, bt_params_key);

    shared_ptr<TLweKey> key = tlwe_keygen(trlweParams_key);
    debug_key = key;
    cout << "time start encrypt X: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<TRGSWMatrix> X = encrypt_X(plaintext_X, *key, N, L_X + 32 + 5, rho);

    cout << "time start encrypt S: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<TRGSWMatrix> S = encrypt_S(plaintext_S, *key, N, L_S + 32 + 5, rho);

    cout << "time start encrypt W: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<TRLWEVector> W = encrypt_individual_trlwe(plaintext_W, *key, N, L_W, tau_W, rho);

    cout << "time start compute A: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<TRLweMatrix> A = compute_A(*X, *S, *W, NA, NA, rho);
    cout << "time end compute A: " << clock() / double(CLOCKS_PER_SEC) << endl;

    TRLWEVector temp(l, A->data[0][0].params);

    mat_RR resp_target;
    resp_target.SetDims(k + 1, l * N / 2);
    clear(resp_target);

    for (int i = 0; i < k + 1; i++) {
        for (int j = 0; j < l * N / 2; j++) {
            for (int kk = 0; kk < n; kk++) {
                resp_target[i][j] += plaintext_X[kk][i] * plaintext_S[kk][j] * plaintext_W[kk];
                //resp_target[i][j] += plaintext_S[kk][j] * plaintext_W[kk];
            }
        }
    }

    for (int i = 0; i < k + 1; i++) {
        for (int j = 0; j < l; j++) {
            copy(temp.data[j], A->data[i][j]);
        }
        vec_RR resp = decrypt_heaan_packed_trlwe(temp, *key, l * N / 2);

        for (int j = 0; j < l * N / 2; j++) {
            cout << resp_target[i][j] << endl;
            cout << resp[j] << endl;
            //EXPECT_LE(log2Diff(resp_target[i][j], resp[j]), tau_S + tau_W + tau_X - rho + 5 + log(N));
             if (j == 100) break; //TODO remove
        }
         break; //TODO remove
    }

}
