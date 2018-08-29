#include <cstdint>
#include <NTL/mat_RR.h>
#include <memory>
#include <fstream>
#include "TLwe.h"
#include "TRLwe.h"
#include "mainalgo.h"
#include "section2_params.h"

NTL_CLIENT;

NTL::mat_RR debug_S;
NTL::mat_RR debug_X;
NTL::vec_RR debug_W;
shared_ptr<TLweKey> debug_key;


int main() {
    // algo parameters (TODO: share these parameters everywhere)
    using namespace section2_params;

    // read the secret key
    cerr << "deserializing the section2 key" << endl;
    ifstream key_in(section2_key_filename);
    assert_dramatically(key_in, "cannot open " << section2_key_filename);
    shared_ptr<TLweKey> key = deserializeTLweKey(key_in, N);
    key_in.close();

    // deserialize bootstrapping key
    cerr << "deserializing level0 key" << endl;
    vector<int64_t> key_lvl0(n_lvl0);
    read_tlwe_key(lvl0_key_filename.c_str(), key_lvl0.data(), n_lvl0);


    // read the input ciphertexts (from section 1)
    // TODO decide the format and read it
    // read the input trlwe

    // ------ p_lvl4
    // ------
    int64_t dummy;
    ifstream p_stream(p_lvl4_filename);
    if (p_stream) {
        istream_read_binary(p_stream, &dummy, sizeof(int64_t));
        assert_dramatically(dummy == algo_n, "wrong size of p");
        auto p_params = deserializeTRLweParams(p_stream);
        assert_dramatically(int64_t(p_params->N) == N);
        assert_dramatically(int64_t(p_params->fixp_params.torus_limbs) == p_limbs);
        assert_dramatically(int64_t(p_params->fixp_params.level_expo) == p_level);
        assert_dramatically(int64_t(p_params->fixp_params.plaintext_expo) == p_plaintext_expo);
        TRLWEVector p_lvl4(algo_n, *p_params);
        for (int64_t i = 0; i < algo_n; i++) {
            deserializeTRLweContent(p_stream, p_lvl4.data[i]);
        }
        cout << "p_lvl4: " << decrypt_individual_trlwe(p_lvl4, *key, algo_n) << endl;
    } else {
        cout << "p_lvl4: " << "absent!" << endl;
    }
    p_stream.close();
}


int main2() {
    int64_t k = 3;
    //int64_t n = 13;
    int64_t n = 253;
    int64_t l = 5;
    //int64_t N = 128;
    int64_t N = 4096;

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

    return 0;

}
