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
    NTL::vec_RR denominator_plaintext;
    NTL::vec_RR numerator_plaintext;

    denominator_plaintext.SetLength(algo_n);
    numerator_plaintext.SetLength(algo_n);


    // read the input ciphertexts (from section 1)
    int64_t length;
    ifstream numerator_stream(numerator_lvl0_filename);
    if (numerator_stream) {
        istream_read_binary(numerator_stream, &length, sizeof(int64_t));
        assert_dramatically(length == algo_n, "wrong size of numerator");
        auto numerator_params = deserializeTRLweParams(numerator_stream);
        assert_dramatically(int64_t(numerator_params->N) == N);
        //assert_dramatically(int64_t(numerator_params->fixp_params.torus_limbs) == p_limbs);
        assert_dramatically(int64_t(numerator_params->fixp_params.level_expo) == numerator_level);
        assert_dramatically(int64_t(numerator_params->fixp_params.plaintext_expo) == numerator_plaintext_expo);
        TRLWEVector numerator(algo_n, *numerator_params);
        for (int64_t i = 0; i < algo_n; i++) {
            deserializeTRLweContent(numerator_stream, numerator.data[i]);
        }
        cout << "numerator: " << decrypt_individual_trlwe(numerator, *key, algo_n) << endl;
        numerator_plaintext = decrypt_individual_trlwe(numerator, *key, algo_n);
    } else {
        cout << "numerator: " << "absent!" << endl;
    }
    numerator_stream.close();

    // serialize numerator plaintext
    ofstream numerator_plaintext_stream("numerator_plaintext_plaintext.bin");
    ostream_write_binary(numerator_plaintext_stream, &length, sizeof(int64_t));
    for (int64_t i = 0; i < algo_n; i++) {
        ostream_write_binary(numerator_plaintext_stream, &numerator_plaintext[i], sizeof(NTL::RR));
    }
    numerator_plaintext_stream.close();

    ifstream denominator_stream(denominator_lvl0_filename);
    if (denominator_stream) {
        istream_read_binary(denominator_stream, &length, sizeof(int64_t));
        assert_dramatically(length == algo_n, "wrong size of denominator");
        auto denominator_params = deserializeTRLweParams(denominator_stream);
        assert_dramatically(int64_t(denominator_params->N) == N);
        assert_dramatically(int64_t(denominator_params->fixp_params.torus_limbs) == denominator_limbs);
        assert_dramatically(int64_t(denominator_params->fixp_params.level_expo) == denominator_level);
        assert_dramatically(int64_t(denominator_params->fixp_params.plaintext_expo) == denominator_plaintext_expo);
        TRLWEVector denominator(algo_n, *denominator_params);
        for (int64_t i = 0; i < algo_n; i++) {
            deserializeTRLweContent(denominator_stream, denominator.data[i]);
        }
        cout << "denominator: " << decrypt_individual_trlwe(denominator, *key, algo_n) << endl;
        denominator_plaintext = decrypt_individual_trlwe(denominator, *key, algo_n);
    } else {
        cout << "denominator: " << "absent!" << endl;
    }
    denominator_stream.close();



    //serialize denumerator plaintext
    ofstream denominator_plaintext_stream("denomerator_plaintext_plaintext.bin");
    ostream_write_binary(denominator_plaintext_stream, &length, sizeof(int64_t));
    for (int64_t i = 0; i < algo_n; i++) {
        ostream_write_binary(denominator_plaintext_stream, &denominator_plaintext[i], sizeof(NTL::RR));
    }
    denominator_plaintext_stream.close();


    NTL::vec_RR stat;
    stat.SetLength(algo_n);
    for (int64_t i = 0; i < algo_n; i++) {
        stat[i] = numerator_plaintext[i] / denominator_plaintext[i];
    }

    //serialize stat
    ofstream stat_stream("stat.bin");
    ostream_write_binary(stat_stream, &length, sizeof(int64_t));
    for (int64_t i = 0; i < algo_n; i++) {
        ostream_write_binary(stat_stream, &stat[i], sizeof(NTL::RR));
    }
    stat_stream.close();

}


