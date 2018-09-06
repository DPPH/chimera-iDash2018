#include <cstdint>
#include <NTL/mat_RR.h>
#include <memory>
#include <fstream>
#include "TLwe.h"
#include "TRLwe.h"
#include "mainalgo.h"
#include "section2_params.h"

NTL_CLIENT;

int main() {
    // algo parameters (TODO: share these parameters everywhere)
    using namespace section2_params;

    // read the secret key
    cerr << "deserializing the section2 key" << endl;
    ifstream key_in(section2_key_filename);
    assert_dramatically(key_in, "cannot open " << section2_key_filename);
    shared_ptr<TLweKey> key = deserializeTLweKey(key_in, N);
    key_in.close();


    NTL::vec_RR denominator_plaintext;
    NTL::vec_RR numerator_plaintext;

    denominator_plaintext.SetLength(algo_m);
    numerator_plaintext.SetLength(algo_m);


    // read the input ciphertexts 
    int64_t length;
    ifstream numerator_stream(numerator_lvl0_filename);
    if (numerator_stream) {
        istream_read_binary(numerator_stream, &length, sizeof(int64_t));
        auto numerator_params = deserializeTRLweParams(numerator_stream);
        store_forever(numerator_params);
        assert_dramatically(int64_t(numerator_params->N) == N);
        TRLWEVector numerator(length, *numerator_params);
        for (int64_t i = 0; i < length; i++) {
            deserializeTRLweContent(numerator_stream, numerator.data[i]);
        }
        numerator_plaintext = decrypt_heaan_packed_trlwe(numerator, *key, algo_m);
        cout << "numerator: " << numerator_plaintext << endl;
       
    } else {
        cout << "numerator: " << "absent!" << endl;
    }
    numerator_stream.close();


    ifstream denominator_stream(denominator_lvl0_filename);
    if (denominator_stream) {
        istream_read_binary(denominator_stream, &length, sizeof(int64_t));
        auto denominator_params = deserializeTRLweParams(denominator_stream);
        store_forever(denominator_params);
        assert_dramatically(int64_t(denominator_params->N) == N);
        TRLWEVector denominator(length, *denominator_params);
        for (int64_t i = 0; i < length; i++) {
            deserializeTRLweContent(denominator_stream, denominator.data[i]);
        }
        denominator_plaintext = decrypt_heaan_packed_trlwe(denominator, *key, algo_m);
        cout << "denominator: " << denominator_plaintext << endl;

    } else {
        cout << "denominator: " << "absent!" << endl;
    }
    denominator_stream.close();


    NTL::vec_RR stat;
    stat.SetLength(algo_m);
    for (int64_t i = 0; i < algo_m; i++) {
        stat[i] = abs(numerator_plaintext[i]) / sqrt(abs(denominator_plaintext[i]));
    }

    //serialize stat
    ofstream stat_stream("stat.txt");
    for (int64_t i = 0; i < algo_m; i++) {
        stat_stream << stat[i] << endl;
    }
    stat_stream.close();

}


