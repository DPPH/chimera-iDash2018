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

//sigmoid function
static RR sigmoid(double x) {
    return to_RR(1) / (to_RR(1) + exp(-to_RR(x)));
}

int main() {
    using namespace section2_params;

    shared_ptr<TRGSWParams> bk_trgsw_params;
    // blind-rotate the input ciphertext

    // read the bk key
    // deserialize bootstrapping key
    cerr << "deserializing bk" << endl;
    ifstream bk_key_in(section1_2_bk_filename);
    int64_t dummy;
    istream_read_binary(bk_key_in, &dummy, sizeof(int64_t));
    assert_dramatically(dummy == n_lvl0);
    bk_trgsw_params = deserializeTRGSWParams(bk_key_in);
    store_forever(bk_trgsw_params);
    TRGSW *bk = new_TRGSW_array(n_lvl0, *bk_trgsw_params);
    for (int j = 0; j < n_lvl0; ++j) {
        deserializeTRGSWContent(bk_key_in, bk[j]);
    }
    bk_key_in.close();

    // read the ks key
    cerr << "deserializing ks" << endl;
    ifstream ks_key_in(section1_2_ks_filename);
    shared_ptr<pubKsKey32> ks_key = deserializepubKsKey32(ks_key_in);
    ks_key_in.close();

    // ------ test vector params
    // ------
    const int64_t N = bk_trgsw_params->N;
    const int64_t Ns2 = N / 2;
    BigTorusParams test_vector_bt_params(2, test_vector_plaintext_expo, test_vector_level);
    TRLweParams test_vector_params(N, test_vector_bt_params);

    // create the test vector corresponding to the sigmoid function
    assert_dramatically(lvl0_ciphertext_modulus == 2 * N);
    vec_RR plaintext_test_vector(INIT_SIZE, N);
    for (int64_t i = 0; i < Ns2; i++) {
        // mapping: -modulus/4 (= -N/2)  corresponds to -2^tau
        // mapping: -modulus/4 (= +N/2)  corresponds to +2^tau
        // modulus: 2N
        plaintext_test_vector[i] = sigmoid(i / double(Ns2) * pow(2., test_vector_plaintext_expo));
    }
    for (int64_t i = 1; i <= Ns2; i++) {
        //tv[N-i]=-tv[-i]
        plaintext_test_vector[N - i] = -sigmoid(-i / double(Ns2) * pow(2., test_vector_plaintext_expo));
    }
    TRLwe sigmoid_test_vector(test_vector_params);
    fixp_trivial(sigmoid_test_vector, plaintext_test_vector, section2_params::default_plaintext_precision);

    // read the input ciphertexts (from section 1)
    // read the input trlwe
    cerr << "reading section 0 ciphertext" << endl;
    int64_t *in_coefs_raw = new int64_t[(n_lvl0 + 1) * algo_n];
    int64_t **in_coefs = new int64_t *[algo_n];
    for (int64_t i = 0; i < algo_n; i++) {
        in_coefs[i] = in_coefs_raw + (n_lvl0 + 1) * i;
    }
    read_tlwe_samples(lvl0_ciphertext_filename.c_str(), in_coefs, algo_n, n_lvl0, 2 * N);

    // ------ p params
    // ------
    BigTorusParams p_bt_params(p_limbs, p_plaintext_expo, p_level);
    TLweParams p_tlwe_params(N, p_bt_params);  // extract after bootstrap
    TRLweParams p_trlwe_params(N, p_bt_params); // param after pubKS

    TRLWEVector p_lvl4(algo_n, p_trlwe_params);

#pragma omp parallel for ordered schedule(dynamic, 1)
    for (int64_t i = 0; i < algo_n; i++) {
#pragma omp critical
        cerr << "bootstrapping p_" << i << endl;
        TRLwe rotated_test_vector(p_trlwe_params);
        TLwe extracted_sigmoid(p_trlwe_params);
        //blind rotate it
        copy(rotated_test_vector, sigmoid_test_vector);
        blind_rotate(rotated_test_vector,
                     in_coefs[i][n_lvl0], in_coefs[i],
                     bk, n_lvl0, p_alpha_bits);
        //extract the constant term
#pragma omp critical
        cerr << "extract p_" << i << endl;
        copy(extracted_sigmoid.getAT(N), rotated_test_vector.a[1].getAT(0)); //constant term
        copy(extracted_sigmoid.getAT(0), rotated_test_vector.a[0].getAT(0));
        for (int64_t j = 1; j < N; j++) {
            neg(extracted_sigmoid.getAT(j), rotated_test_vector.a[0].getAT(N - j));
        }
        //pubks it to p_lvl4
#pragma omp critical
        cerr << "pubKS p_" << i << endl;
        pubKS32(p_lvl4.data[i], extracted_sigmoid, *ks_key, p_limbs);
    }

    // free resources for ks_key and bk_key
    ks_key = nullptr;
    delete_TRGSW_array(n_lvl0, bk);
    delete[] in_coefs;
    delete[] in_coefs_raw;

    // serialize p (lvl 4)
    ofstream p_stream("p_lvl4.bin");
    ostream_write_binary(p_stream, &p_lvl4.length, sizeof(int64_t));
    serializeTRLweParams(p_stream, p_lvl4.data[0].params);
    for (int64_t i = 0; i < algo_n; i++) {
        serializeTRLweContent(p_stream, p_lvl4.data[i]);
    }
    p_stream.close();

    // ------------------
    // compute W (lvl 3)

    // -----------------
    // compute numerator (lvl 0)   (requires enc. y of S)

    // ------------------
    // compute A (lvl 1)           (requires enc. of S and X)


    //free S and X here

    // -------------------
    // denom_1 proportional to row 0 of A

    // -------------------
    // denom 2 = 4*norms2(A)
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
