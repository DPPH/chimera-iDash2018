#include <cstdint>
#include <NTL/mat_RR.h>
#include <memory>
#include <fstream>
#include <cstring>
#include "TLwe.h"
#include "TRLwe.h"
#include "mainalgo.h"
#include "section2_params_temporal.h"

NTL_CLIENT;

int main() {
    // algo parameters (TODO: share these parameters everywhere)
    using namespace section2_params_temporal;

    // read the secret key
    cerr << "deserializing the section2 key" << endl;
    ifstream key_in(section2_key_filename);
    assert_dramatically(key_in, "cannot open " << section2_key_filename);
    shared_ptr<TLweKey> key = deserializeTLweKey(key_in, N);
    key_in.close();


    int length;
    //deserialize normally
    ifstream numerator_stream("numerator_lvl0.bin");
        istream_read_binary(numerator_stream, &length, sizeof(int64_t));
    shared_ptr<TRLweParams> params = deserializeTRLweParams(numerator_stream);
    store_forever(params);
    TRLWEVector numeratorTmp(length, *params);
    for (int64_t i = 0; i < numeratorTmp.length; i++) {
        deserializeTRLweContent(numerator_stream, numeratorTmp.data[i]);
        }
    numerator_stream.close();

    TRLWEVector *numerator = &numeratorTmp;


    cerr << "Re-serialize numerator for section 3" << endl;
    //numerator is already at level 1, so no need to shift
    ofstream N0_stream("num_sec3.bin");
    int64_t tmp = N;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //key size
    tmp = numerator->length;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //number of ciphertexts
    tmp = algo_m;
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t)); //unpacked length
    tmp = numerator->data[0].params.fixp_params.plaintext_expo + 1; //scaling factor
    ostream_write_binary(N0_stream, &tmp, sizeof(int64_t));
    for (int64_t i = 0; i < numerator->length; i++) {
        for (int64_t k = 0; k < N; k++) {
            ostream_write_binary(N0_stream, ((char *) numerator->data[i].a[0].getAT(k).limbs_end) - 4, 4);
        }
        for (int64_t k = 0; k < N; k++) {
            ostream_write_binary(N0_stream, ((char *) numerator->data[i].a[1].getAT(k).limbs_end) - 4, 4);
        }
    }
    N0_stream.close();




    //deserialize sec3
    int64_t sec3_N;
    int64_t sec3_length;
    int64_t sec3_algo_m;
    int64_t sec3_scaling;
    ifstream N0_stream1("num_sec3.bin");
    istream_read_binary(N0_stream1, &sec3_N, sizeof(int64_t)); //key size
    assert_dramatically(sec3_N == N);
    istream_read_binary(N0_stream1, &sec3_length, sizeof(int64_t)); //number of ciphertexts
    istream_read_binary(N0_stream1, &sec3_algo_m, sizeof(int64_t)); //unpacked length
    istream_read_binary(N0_stream1, &sec3_scaling, sizeof(int64_t));

    int32_t *sec3_a = new int32_t[sec3_N];
    int32_t *sec3_b = new int32_t[sec3_N];
    BigTorusParams test_bt_params(1, sec3_scaling - 1, 1);
    TRLweParams test_params(sec3_N, test_bt_params);
    //numerator->data[0].params.fixp_params.plaintext_expo + 1; //scaling factor
    for (int64_t i = 0; i < sec3_length; i++) {
        istream_read_binary(N0_stream1, sec3_a, sec3_N * sizeof(int32_t));
        istream_read_binary(N0_stream1, sec3_b, sec3_N * sizeof(int32_t));
        //decrypt it
        TRLwe test(test_params);
        zero(test);
        for (int64_t k = 0; k < sec3_N; k++) {
            to_torus(test.a[0].getAT(k), to_RR(double(sec3_a[k])) / power2_RR(32));
            //memcpy(test.a[0].getAT(k).limbs_end-4, sec3_a+k, 4);
        }
        for (int64_t k = 0; k < sec3_N; k++) {
            to_torus(test.a[1].getAT(k), to_RR(double(sec3_b[k])) / power2_RR(32));
            //memcpy(test.a[1].getAT(k).limbs_end-4, sec3_b+k, 4);
        }
        vec_RR reps = fixp_decrypt(test, *key);
        cout << reps << endl;
    }
    N0_stream1.close();

    // deserialize A (lvl1)
    ifstream Ain("A_lvl1.bin");
    int64_t Arows;
    int64_t Acols;
    istream_read_binary(Ain, &Arows, sizeof(int64_t));
    istream_read_binary(Ain, &Acols, sizeof(int64_t));
    shared_ptr<TRLweParams> Aparams = deserializeTRLweParams(Ain);
    TRLweMatrix AA(Arows, Acols, *Aparams);
    TRLweMatrix *A = &AA;

    for (int64_t i = 0; i < A->rows; i++) {
        for (int64_t j = 0; j < A->cols; j++) {
            deserializeTRLweContent(Ain, A->data[i][j]);
        }
    }
    Ain.close();


    // serialize A (lvl 0) -> section 3
    cerr << "Re-Serialize A for section 3" << endl;
    BigTorusParams a_lvl0_bt_params(1, 0, 1);
    TRLweParams a_lvl0_trlwe_params(N, a_lvl0_bt_params);
    int64_t rsShift = A->data[0][0].params.fixp_params.level_expo - 1;
    TRLwe tmpA(a_lvl0_trlwe_params);

    ofstream A0_stream("A_sec3.bin");
    //int64_t tmp;
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
            lshift(tmpA, A->data[i][j], rsShift);
            for (int64_t k = 0; k < N; k++) {
                ostream_write_binary(A0_stream, ((char *) tmpA.a[0].getAT(k).limbs_end) - 4, 4);
            }
            for (int64_t k = 0; k < N; k++) {
                ostream_write_binary(A0_stream, ((char *) tmpA.a[1].getAT(k).limbs_end) - 4, 4);
            }
        }
    }
    A0_stream.close();

    //Re-deserialize to test
    ifstream A0_stream2("A_sec3.bin");
    //int64_t tmp;
    istream_read_binary(A0_stream2, &tmp, sizeof(int64_t)); //key size
    assert_dramatically(tmp == N);
    istream_read_binary(A0_stream2, &tmp, sizeof(int64_t)); //rows
    assert_dramatically(tmp == A->rows);
    istream_read_binary(A0_stream2, &tmp, sizeof(int64_t)); //rows
    assert_dramatically(tmp == A->cols);
    istream_read_binary(A0_stream2, &tmp, sizeof(int64_t)); //rows
    assert_dramatically(tmp == algo_m);
    istream_read_binary(A0_stream2, &tmp, sizeof(int64_t)); //rows
    assert_dramatically(tmp == A->data[0][0].params.fixp_params.plaintext_expo + 1);
    BigTorusParams AtestBtParams2(1, A->data[0][0].params.fixp_params.plaintext_expo, 1);
    TRLweParams AtestParams2(N, AtestBtParams2);
    for (int64_t i = 0; i < A->rows; i++) {
        for (int64_t j = 0; j < A->cols; j++) {
            TRLwe test(AtestParams2);
            istream_read_binary(A0_stream2, sec3_a, sec3_N * sizeof(int32_t));
            istream_read_binary(A0_stream2, sec3_b, sec3_N * sizeof(int32_t));
            zero(test);
            for (int64_t k = 0; k < N; k++) {
                to_torus(test.a[0].getAT(k), to_RR(double(sec3_a[k])) / power2_RR(32));
                //memcpy(test.a[0].getAT(k).limbs_end-4, sec3_a+k, 4);
            }
            for (int64_t k = 0; k < N; k++) {
                to_torus(test.a[1].getAT(k), to_RR(double(sec3_b[k])) / power2_RR(32));
                //memcpy(test.a[1].getAT(k).limbs_end-4, sec3_b+k, 4);
            }
            vec_RR reps = fixp_decrypt(test, *key);
            cout << "A" << i << j << " " << reps << endl;
        }
    }
    A0_stream2.close();

}


