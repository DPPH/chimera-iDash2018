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

#define DEBUG_MODE

int main() {
    using namespace section2_params_temporal;


    cerr << "Deserializing the section2 key (N = " << N << ")" << endl;
    ifstream key_in(section2_key_filename);
    assert_dramatically(key_in, "cannot open " << section2_key_filename);
    shared_ptr<TLweKey> key = deserializeTLweKey(key_in, N);
    assert_dramatically(key_in, "problem in key deserialization");
    key_in.close();

    mat_RR X(INIT_SIZE, algo_n, algo_k + 1);
    mat_RR S(INIT_SIZE, algo_n, algo_m);
    vec_RR y(INIT_SIZE, algo_n);

    cerr << "reading X,y..." << endl;
    fill_matrix_Xy(X, y);
    cerr << "reading S..." << endl;
    fill_matrix_S(S);


    cerr << "encrypt y..." << endl;
    shared_ptr<TRLWEVector> ciphertext_y = encrypt_individual_trlwe(y, *key, N, y_level, y_plaintext_expo,
                                                                    section2_params_temporal::default_plaintext_precision);
    cerr << "encrypt X..." << endl;
    shared_ptr<TRGSWMatrix> ciphertext_X = encrypt_X(X, *key, N, X_alpha_bits,
                                                     section2_params_temporal::default_plaintext_precision);

    cerr << "encrypt S..." << endl;
    shared_ptr<TRGSWMatrix> ciphertext_S = encrypt_S_temporal(S, *key, N, S_alpha_bits, 1); // S is binary

#ifdef DEBUG_MODE
    cerr << "y :" << y << endl;
    cerr << "encrypted_y :" << decrypt_individual_trlwe(*ciphertext_y, *key, algo_n) << endl;
#endif


    // serialize y (lvl 2)
    ofstream y_stream(y_lvl2_filename);
    ostream_write_binary(y_stream, &ciphertext_y->length, sizeof(int64_t));
    serializeTRLweParams(y_stream, ciphertext_y->data[0].params);
    for (int64_t i = 0; i < algo_n; i++) {
        serializeTRLweContent(y_stream, ciphertext_y->data[i]);
    }
    y_stream.close();

    // serialize S (lvl 3)
    ofstream S_stream(S_lvl3_filename);
    ostream_write_binary(S_stream, &ciphertext_S->rows, sizeof(int64_t));
    ostream_write_binary(S_stream, &ciphertext_S->cols, sizeof(int64_t));
    serializeTRGSWParams(S_stream, ciphertext_S->data[0][0].params);

    for (int64_t i = 0; i < ciphertext_S->rows; i++) {
        for (int64_t j = 0; j < ciphertext_S->cols; j++) {
            serializeTRGSWContent(S_stream, ciphertext_S->data[i][j]);
        }
    }
    S_stream.close();

    // serialize X (lvl 2)
    ofstream X_stream(X_lvl2_filename);
    ostream_write_binary(X_stream, &ciphertext_X->rows, sizeof(int64_t));
    ostream_write_binary(X_stream, &ciphertext_X->cols, sizeof(int64_t));
    serializeTRGSWParams(X_stream, ciphertext_X->data[0][0].params);

    for (int64_t i = 0; i < ciphertext_X->rows; i++) {
        for (int64_t j = 0; j < ciphertext_X->cols; j++) {
            serializeTRGSWContent(X_stream, ciphertext_X->data[i][j]);
        }
    }
    X_stream.close();


}

