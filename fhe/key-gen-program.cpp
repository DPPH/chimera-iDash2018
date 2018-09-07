#include <cstdint>
#include <memory>
#include <fstream>
#include <iostream>
#include <cassert>
#include "BigTorus.h"
#include "TRGSW.h"
#include "mainalgo.h"
#include "section2_params.h"

using namespace std;

int main() {
    using namespace section2_params;

    // Generation of Bootstrapping key (s)

    BigTorusParams bk_bt_params(bk_nblimbs);
    TRGSWParams bk_trgswParams(N, bk_bt_params);

    BigTorusPolynomial phase(N, bk_bt_params);

    TRLwe reps(bk_trgswParams);
    TRGSW *bk = new_TRGSW_array(n_lvl0, bk_trgswParams);

    shared_ptr<TLweKey> key = tlwe_keygen(bk_trgswParams);

    vector<int64_t> s(n_lvl0);
    read_tlwe_key(lvl0_key_filename.c_str(), s.data(), n_lvl0);

    cout << "start encrypt bootstrapping key: " << clock() / double(CLOCKS_PER_SEC) << endl;

    int k = 1;
#pragma omp parallel for

    for (int i = 0; i < n_lvl0; i++) {

        int_encrypt(bk[i], s[i], *key, bk_alpha_bits + TRGSWParams::Bgbits);
#pragma omp critical
        {
            printf("%3d/%3ld\r", k++, long(n_lvl0));
            fflush(stdout);
        }
    }
    printf("\n");

    cout << "end encrypt bootstrapping key: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //Generation KS key (ks_key)

    int64_t ks_N_in = bk_trgswParams.N;
    int64_t ks_N_out = N;

    BigTorusParams &ks_bt_params_in = bk_bt_params;
    BigTorusParams ks_bt_params_out(ks_nblimbs_out);
    TLweParams tlwe_params_in(ks_N_in, ks_bt_params_in);
    TRLweParams trlwe_params_out(ks_N_out, ks_bt_params_out);

    shared_ptr<TLweKey> key_in = key;
    shared_ptr<TLweKey> key_out = key;

    cout << "start keygen key switch key at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<pubKsKey32> ks_key = ks_keygen32(
            trlwe_params_out, tlwe_params_in,
            *key_in, *key_out, ks_out_alpha_bits);
    cout << "end keygen key switch key at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //Generation RK key (rk)

    BigTorusParams bt_params_rk(rk_nblimbs);
    TRGSWParams trgswParams_rk(N, bt_params_rk);
    TRGSW rk(trgswParams_rk);
    cout << "start keygen rk key at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    intPoly_encrypt(rk, key->key, *key, rk_alpha_bits + TRGSWParams::Bgbits);
    cout << "end keygen rk key at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //serialize bootstrapping key
    ofstream bk_key_out(section1_2_bk_filename);
    ostream_write_binary(bk_key_out, &n_lvl0, sizeof(int64_t));
    serializeTRGSWParams(bk_key_out, bk_trgswParams);
    for (int j = 0; j < n_lvl0; ++j) {
        serializeTRGSWContent(bk_key_out, bk[j]);
    }
    bk_key_out.close();

    //serialize the secret key
    ofstream secret_key_out(section2_key_filename);
    serializeTLweKey(secret_key_out, *key, N);
    secret_key_out.close();

    //serialize key switch key
    ofstream ks_key_out(section1_2_ks_filename);
    serializepubKsKey32(ks_key_out, *ks_key);
    ks_key_out.close();

    //serialize relian. key
    ofstream rk_key_out(section2_rk_filename);
    serializeTRGSW(rk_key_out, rk);
    rk_key_out.close();

}
