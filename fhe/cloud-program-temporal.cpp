#include <cstdint>
#include <NTL/mat_RR.h>
#include <memory>
#include <fstream>
#include "TLwe.h"
#include "TRLwe.h"
#include "mainalgo.h"
#include "section2_params_temporal.h"
#include "arithmetic.h"

NTL_CLIENT;

//sigmoid function
RR sigmoid(double x) {
    return to_RR(1) / (to_RR(1) + exp(-to_RR(x)));
}

int64_t compute_public_exponent(const vec_RR &v) {
    double pubExpo = -INFINITY;
    for (int64_t i = 0; i < v.length(); i++) {
        if (v[i] == 0) continue;
        pubExpo = std::max(pubExpo, to_double(log(abs(v[i]))) / log(2));
    }
    cerr << "computed plaintext expo: " << pubExpo << endl;
    return int64_t(ceil(pubExpo));
}

int64_t compute_public_exponent(const mat_RR &A) {
    double pubExpo = -INFINITY;
    for (int64_t i = 0; i < A.NumRows(); i++) {
        for (int64_t j = 0; j < A.NumCols(); j++) {
            if (A[i][j] == 0) continue;
            pubExpo = std::max(pubExpo, to_double(log(abs(A[i][j]))) / log(2.));
        }
    }
    cerr << "computed plaintext expo: " << pubExpo << endl;
    return int64_t(ceil(pubExpo));
}

void print_difference(const vec_RR &actual, const vec_RR &expected, const string &name) {
    RR max_diff;
    max_diff = 0;
    RR avg_diff;
    avg_diff = 0;
    ofstream ofs(name + string(".dat"));
    for (int64_t i = 0; i < actual.length(); i++) {
        ofs << i << " " << actual[i] << " " << expected[i] << endl;
        RR diff = abs(actual[i] - expected[i]);
        if (diff > max_diff) max_diff = diff;
        avg_diff += diff;
    }
    ofs.close();
    avg_diff /= actual.length();
    cerr << "Distance " << name << " max-error: " << max_diff << " avg-error: " << avg_diff << endl;
}

void print_difference(const mat_RR &actual, const mat_RR &expected, const string &name) {
    RR max_diff;
    max_diff = 0;
    RR avg_diff;
    avg_diff = 0;
    ofstream ofs(name + string(".dat"));
    for (int64_t i = 0; i < actual.NumRows(); i++) {
        for (int64_t j = 0; j < actual.NumCols(); j++) {
            ofs << i << " " << j << " " << actual[i][j] << " " << expected[i][j] << endl;
            RR diff = abs(actual[i][j] - expected[i][j]);
            if (diff > max_diff) max_diff = diff;
            avg_diff += diff;
        }
        ofs << endl;
    }
    ofs.close();
    avg_diff /= (actual.NumRows() * actual.NumCols());
    cerr << "Distance " << name << " max-error: " << max_diff << " avg-error: " << avg_diff << endl;
}

#define DEBUG_MODE

int main() {
    using namespace section2_params_temporal;

#ifdef DEBUG_MODE
    cerr << "DEBUG!!! reading lvl0 key...";
    vector<int64_t> s(n_lvl0);
    read_tlwe_key(lvl0_key_filename.c_str(), s.data(), n_lvl0);

    cerr << "DEBUG!!! deserializing the section2 key" << endl;
    ifstream key_in(section2_key_filename);
    assert_dramatically(key_in, "cannot open " << section2_key_filename);
    shared_ptr<TLweKey> key = deserializeTLweKey(key_in, N);
    assert_dramatically(key_in, "problem in keyswitch key");
    key_in.close();

    mat_RR plain_S(INIT_SIZE, algo_n, algo_m);
    mat_RR plain_X(INIT_SIZE, algo_n, algo_k + 1);
    vec_RR plain_y(INIT_SIZE, algo_n);
    fill_matrix_S(plain_S);
    fill_matrix_Xy(plain_X, plain_y);
#endif


//#define FAKE_BOOTSTRAPPING
#ifdef FAKE_BOOTSTRAPPING
    // ------ test vector and p params
    // ------
    cerr << "DEBUG: faking the bootstrapping completely!!" << endl;
    BigTorusParams p_bt_params(p_limbs, p_plaintext_expo, p_level);
    TLweParams p_tlwe_params(N, p_bt_params);  // extract after bootstrap
    TRLweParams p_trlwe_params(N, p_bt_params); // param after pubKS
    vec_RR pp(INIT_SIZE, algo_n);
    for (int64_t i = 0; i < algo_n; i++) {
        pp[i] = random_RR();
    }
    shared_ptr<TRLWEVector> p_lvl4p = encrypt_individual_trlwe(pp, *key, N, p_level, p_plaintext_expo,
                                                               section2_params_temporal::default_plaintext_precision);
    TRLWEVector &p_lvl4 = *p_lvl4p;
#else
    // read the bk key
    // deserialize bootstrapping key
    cerr << "deserializing bk" << endl;
    shared_ptr<TRGSWParams> bk_trgsw_params;
    ifstream bk_key_in(section1_2_bk_filename);
    int64_t dummy;
    istream_read_binary(bk_key_in, &dummy, sizeof(int64_t));
    assert_dramatically(dummy == n_lvl0);
    bk_trgsw_params = deserializeTRGSWParams(bk_key_in);
    store_forever(bk_trgsw_params);
    TRGSW *bk = new_TRGSW_array(n_lvl0, *bk_trgsw_params);
    for (int j = 0; j < n_lvl0; ++j) {
        deserializeTRGSWContent(bk_key_in, bk[j]);
        assert_dramatically(bk_key_in, "problem in the serialization file");
    }
    bk_key_in.close();

    // read the ks key
    cerr << "deserializing ks" << endl;
    ifstream ks_key_in(section1_2_ks_filename);
    shared_ptr<pubKsKey32> ks_key = deserializepubKsKey32(ks_key_in);
    ks_key_in.close();

    // ------ test vector and p params
    // ------
    BigTorusParams p_bt_params(p_limbs, p_plaintext_expo, p_level);
    TLweParams p_tlwe_params(N, p_bt_params);  // extract after bootstrap
    TRLweParams p_trlwe_params(N, p_bt_params); // param after pubKS

    const int64_t Ns2 = N / 2;

    // create the test vector corresponding to the sigmoid function
    assert_dramatically(lvl0_ciphertext_modulus == 2 * N);
    vec_RR plaintext_test_vector(INIT_SIZE, N);
    for (int64_t i = 0; i < Ns2; i++) {
        // mapping: -modulus/4 (= -N/2)  corresponds to -2^tau
        // mapping: -modulus/4 (= +N/2)  corresponds to +2^tau
        // modulus: 2N
        plaintext_test_vector[i] = sigmoid(i / double(Ns2) * pow(2., lvl0_ciphertext_plaintext_expo)) - 0.5;
    }
    for (int64_t i = 1; i <= Ns2; i++) {
        //tv[N-i]=-tv[-i]
        plaintext_test_vector[N - i] = -sigmoid(-i / double(Ns2) * pow(2., lvl0_ciphertext_plaintext_expo)) + 0.5;
    }
    BigTorus sigmoid_offset(p_bt_params);
    to_fixP(sigmoid_offset, to_RR(0.5));
    TRLwe sigmoid_test_vector(p_trlwe_params);
    fixp_trivial(sigmoid_test_vector, plaintext_test_vector, section2_params_temporal::default_plaintext_precision);

    // read the input ciphertexts (from section 1)
    // read the input trlwe
    cerr << "reading section 0 ciphertext" << endl;
    int64_t *in_coefs_raw = new int64_t[(n_lvl0 + 1) * algo_n];
    int64_t **in_coefs = new int64_t *[algo_n];
    for (int64_t i = 0; i < algo_n; i++) {
        in_coefs[i] = in_coefs_raw + (n_lvl0 + 1) * i;
    }
    read_tlwe_samples(lvl0_ciphertext_filename.c_str(), in_coefs, algo_n, n_lvl0, 2 * N);

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
                     bk, n_lvl0, p_alpha_bits + 10);
        //extract the constant term
#pragma omp critical
        cerr << "extract p_" << i << endl;
        add(extracted_sigmoid.getAT(N), rotated_test_vector.a[1].getAT(0), sigmoid_offset); //constant term
        copy(extracted_sigmoid.getAT(0), rotated_test_vector.a[0].getAT(0));
        for (int64_t j = 1; j < N; j++) {
            neg(extracted_sigmoid.getAT(j), rotated_test_vector.a[0].getAT(N - j));
        }

        //pubks it to p_lvl4
#pragma omp critical
        cerr << "pubKS p_" << i << endl;
        pubKS32(p_lvl4.data[i], extracted_sigmoid, *ks_key, p_limbs);

#ifdef DEBUG_MODE
        //input phase
        int64_t in_phase = in_coefs[i][n_lvl0];
        for (int j = 0; j < n_lvl0; j++) in_phase -= in_coefs[i][j] * s[j];
        in_phase = ((in_phase % lvl0_ciphertext_modulus) + lvl0_ciphertext_modulus) % lvl0_ciphertext_modulus;
        //expected output phase
        RR expected_phase;
        vec_RR actual_test_vector = fixp_decrypt(sigmoid_test_vector, *key);
        RR actual_phase_trivial;
        if (in_phase < N) {
            expected_phase = plaintext_test_vector[in_phase];
            actual_phase_trivial = actual_test_vector[in_phase];
        } else {
            expected_phase = -plaintext_test_vector[in_phase - N];
            actual_phase_trivial = -actual_test_vector[in_phase - N];
        }
        //actual phase
        RR actual_phase_blindrotate = fixp_decrypt_number(rotated_test_vector, *key);
        RR actual_phase_extract = slot_decrypt(extracted_sigmoid, *key);
        RR actual_phase_ks = fixp_decrypt_number(p_lvl4.data[i], *key);
#pragma omp critical
        {
            cout << "at index i............: " << i << endl;
            cout << "in phase..............: " << in_phase << endl;
            cout << "expected output phase.: " << expected_phase << endl;
            cout << "actual phase trivial..: " << actual_phase_trivial << endl;
            cout << "actual blindrotate....: " << actual_phase_blindrotate << endl;
            cout << "actual extract........: " << actual_phase_extract << endl;
            cout << "actual ks.............: " << actual_phase_ks << endl;
        }
#endif //DEBUG_MODE
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
#endif //FAKE_BOOTSTRAPPING

#ifdef DEBUG_MODE
    vec_RR decrypted_p = decrypt_individual_trlwe(p_lvl4, *key, algo_n);
    cerr << "DEBUG decrypt p: " << decrypted_p << endl;
    print_difference(decrypted_p, decrypted_p, "p");
    assert_weakly(compute_public_exponent(decrypted_p) <= p_lvl4.data[0].params.fixp_params.plaintext_expo);
#endif //DEBUG_MODE

    // ------------------
    // compute W (lvl 3) w = p(1-p)
    cerr << "deserializing rk" << endl;
    ifstream rk_key_in(section2_rk_filename);
    shared_ptr<TRGSW> rk = deserializeTRGSW(rk_key_in);
    rk_key_in.close();
    shared_ptr<TRLWEVector> w_lvl3 = compute_w(p_lvl4, *rk, section2_params_temporal::w_level,
                                               section2_params_temporal::w_plaintext_expo,
                                               section2_params_temporal::default_plaintext_precision);

    // serialize w (lvl 3)
    ofstream w_stream("w_lvl3.bin");
    ostream_write_binary(w_stream, &w_lvl3->length, sizeof(int64_t));
    serializeTRLweParams(w_stream, w_lvl3->data[0].params);
    for (int64_t i = 0; i < algo_n; i++) {
        serializeTRLweContent(w_stream, w_lvl3->data[i]);
    }
    w_stream.close();

#ifdef DEBUG_MODE
    vec_RR expected_W(INIT_SIZE, algo_n);
    vec_RR decrypted_W = decrypt_individual_trlwe(*w_lvl3, *key, algo_n);
    for (int i = 0; i < algo_n; i++) expected_W[i] = decrypted_p[i] - decrypted_p[i] * decrypted_p[i];
    cerr << "DEBUG decrypt w: " << decrypted_W << endl;
    cerr << "DEBUG expect w: " << expected_W << endl;
    print_difference(decrypted_W, expected_W, "W");
    assert_weakly(compute_public_exponent(decrypted_W) <= w_lvl3->data[0].params.fixp_params.plaintext_expo);
#endif //DEBUG_MODE



    // -----------------
    // compute numerator (lvl 0)   (requires enc. y of S)
    // read y
    cerr << "deserializing y" << endl;
    ifstream y_in(y_lvl2_filename);
    int64_t y_length;
    istream_read_binary(y_in, &y_length, sizeof(int64_t));

    shared_ptr<TRLweParams> y_params = deserializeTRLweParams(y_in);
    TRLWEVector y(y_length, *y_params);

    for (int64_t i = 0; i < algo_n; i++) {
        deserializeTRLweContent(y_in, y.data[i]);
    }
    y_in.close();

#ifdef DEBUG_MODE
    vec_RR decrypted_y = decrypt_individual_trlwe(y, *key, algo_n);
    cerr << "DEBUG decrypt y: " << decrypt_individual_trlwe(y, *key, algo_n) << endl;
    cerr << "DEBUG expect y: " << plain_y << endl;
#endif //DEBUG_MODE

    shared_ptr<TRLWEVector> ymp = substract_ind_TRLWE(y, p_lvl4, y_level, y_plaintext_expo,
                                                      section2_params_temporal::default_plaintext_precision);
#ifdef DEBUG_MODE
    vec_RR expected_ymp = decrypted_y - decrypted_p;
    vec_RR decrypted_ymp = decrypt_individual_trlwe(*ymp, *key, algo_n);
    cerr << "DEBUG decrypt ymp: " << decrypted_ymp << endl;
    cerr << "DEBUG expectd ymp: " << expected_ymp << endl;
    print_difference(decrypted_ymp, expected_ymp, "ymp");
    assert_weakly(compute_public_exponent(decrypted_ymp) <= ymp->data[0].params.fixp_params.plaintext_expo);
#endif //DEBUG_MODE

    // read S
    cerr << "deserializing S" << endl;
    ifstream S_in(S_lvl3_filename);
    int64_t S_rows;
    int64_t S_cols;
    istream_read_binary(S_in, &S_rows, sizeof(int64_t));
    istream_read_binary(S_in, &S_cols, sizeof(int64_t));


    shared_ptr<TRGSWParams> S_params = deserializeTRGSWParams(S_in);
    TRGSWMatrix *S = new TRGSWMatrix(S_rows, S_cols, *S_params);
    for (int64_t i = 0; i < S_rows; i++) {
        for (int64_t j = 0; j < S_cols; j++) {
            deserializeTRGSWContent(S_in, S->data[i][j]);
        }
    }
    S_in.close();
    //(y-p)*S

//#ifdef DEBUG_MODE
//    cerr << "DEBUG decrypt S: " << decrypt_heaan_packed_trlwe(*S, *key, p_lvl4.length) << endl;
//#endif //DEBUG_MODE



    shared_ptr<TRLWEVector> numerator = vec_mat_prod(*ymp, *S, numerator_level, numerator_plaintext_expo,
                                                     section2_params_temporal::default_plaintext_precision);


    // serialize numerator (lvl 0)
    ofstream numerator_stream("numerator_lvl0.bin");
    ostream_write_binary(numerator_stream, &numerator->length, sizeof(int64_t));
    serializeTRLweParams(numerator_stream, numerator->data[0].params);
    for (int64_t i = 0; i < numerator->length; i++) {
        serializeTRLweContent(numerator_stream, numerator->data[i]);
    }
    numerator_stream.close();

#ifdef DEBUG_MODE
    vec_RR expected_numerator = decrypted_ymp * plain_S;
    vec_RR decrypted_numerator = decrypt_temporal_packed_trlwe(*numerator, *key, algo_m);
    cerr << "DEBUG decrypt numerator: " << decrypted_numerator << endl;
    cerr << "DEBUG expectd numerator: " << expected_numerator << endl;
    print_difference(decrypted_numerator, expected_numerator, "numerator");
    assert_weakly(compute_public_exponent(decrypted_numerator) <= numerator->data[0].params.fixp_params.plaintext_expo);
#endif //DEBUG_MODE


    // ------------------
    // compute A (lvl 1)           (requires enc. of S and X)

    // read X
    cerr << "deserializing X" << endl;
    ifstream X_in(X_lvl2_filename);
    int64_t X_rows;
    int64_t X_cols;
    istream_read_binary(X_in, &X_rows, sizeof(int64_t));
    istream_read_binary(X_in, &X_cols, sizeof(int64_t));

    shared_ptr<TRGSWParams> X_params = deserializeTRGSWParams(X_in);
    TRGSWMatrix *X = new TRGSWMatrix(X_rows, X_cols, *X_params);
    for (int64_t i = 0; i < X_rows; i++) {
        for (int64_t j = 0; j < X_cols; j++) {
            deserializeTRGSWContent(X_in, X->data[i][j]);
        }
    }
    X_in.close();

    shared_ptr<TRLweMatrix> A = compute_A(*X, *S, *w_lvl3, A_level, A_plaintext_expo,
                                          section2_params_temporal::default_plaintext_precision);

    // serialize A (lvl 1)
    ofstream A_stream("A_lvl1.bin");
    ostream_write_binary(numerator_stream, &numerator->length, sizeof(int64_t));

    ostream_write_binary(A_stream, &A->rows, sizeof(int64_t));
    ostream_write_binary(A_stream, &A->cols, sizeof(int64_t));
    serializeTRLweParams(A_stream, A->data[0][0].params);

    for (int64_t i = 0; i < A->rows; i++) {
        for (int64_t j = 0; j < A->cols; j++) {
            serializeTRLweContent(A_stream, A->data[i][j]);
        }
    }

    A_stream.close();

#ifdef DEBUG_MODE
    mat_RR expected_A(INIT_SIZE, algo_k + 1, algo_m);
    clear(expected_A);
    for (int64_t i = 0; i < algo_k + 1; i++) {
        for (int64_t j = 0; j < algo_m; j++) {
            for (int64_t k = 0; k < algo_n; k++) {
                expected_A[i][j] += decrypted_W[k] * plain_X[k][i] * plain_S[k][j];
            }
        }
    }
    mat_RR decrypted_A = decrypt_temporal_packed_trlwe(*A, *key, algo_m);
    cerr << "DEBUG decrypt A: " << decrypted_A << endl;
    cerr << "DEBUG expectd A: " << expected_A << endl;
    print_difference(decrypted_A, expected_A, "A");
    assert_weakly(compute_public_exponent(decrypted_A) <= A->data[0][0].params.fixp_params.plaintext_expo);
#endif //DEBUG_MODE

    //free S and X here
    S = nullptr;
    X = nullptr;

    // serialize A (lvl 0) -> section 3
    cerr << "Serialize A for section 3" << endl;
    BigTorusParams a_lvl0_bt_params(1, 0, 1);
    TRLweParams a_lvl0_trlwe_params(N, a_lvl0_bt_params);
    int64_t rsShift = A->data[0][0].params.fixp_params.level_expo - 1;
    TRLwe tmpA(a_lvl0_trlwe_params);

    ofstream A0_stream("A_sec3.bin");
    int64_t tmp;
    tmp = N;
    ostream_write_binary(A0_stream, &tmp, sizeof(int64_t)); //key size
    tmp = A->rows;
    ostream_write_binary(A0_stream, &tmp, sizeof(int64_t)); //rows
    tmp = A->cols;
    ostream_write_binary(A0_stream, &tmp, sizeof(int64_t)); //cols
    tmp = algo_m;
    ostream_write_binary(A0_stream, &tmp, sizeof(int64_t)); //unpacked cols
    tmp = A->data[0][0].params.fixp_params.plaintext_expo + 1; //scaling factor
    ostream_write_binary(A0_stream, &tmp, sizeof(int64_t));
    for (int64_t i = 0; i < A->rows; i++) {
        for (int64_t j = 0; j < A->cols; j++) {
            lshift(tmpA, A->data[0][0], rsShift);
            for (int64_t k = 0; k < N; k++) {
                ostream_write_binary(A0_stream, tmpA.a[0].getAT(k).limbs_end - 4, 4);
            }
            for (int64_t k = 0; k < N; k++) {
                ostream_write_binary(A0_stream, tmpA.a[1].getAT(k).limbs_end - 4, 4);
            }
        }
    }
    A0_stream.close();

    cerr << "Serialize numerator for section 3" << endl;
    //numerator is already at level 1, so no need to shift
    ofstream N0_stream("num_sec3.bin");
    tmp = N;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //key size
    tmp = numerator->length;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //number of ciphertexts
    tmp = algo_m;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //unpacked length
    tmp = numerator->data[0].params.fixp_params.plaintext_expo + 1; //scaling factor
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t));
    for (int64_t i = 0; i < numerator->length; i++) {
        for (int64_t k = 0; k < N; k++) {
            ostream_write_binary(N0_stream, numerator->data[i].a[0].getAT(k).limbs_end - 4, 4);
        }
        for (int64_t k = 0; k < N; k++) {
            ostream_write_binary(N0_stream, numerator->data[i].a[1].getAT(k).limbs_end - 4, 4);
        }
    }
    N0_stream.close();


}
