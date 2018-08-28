#include <cstdint>
#include <memory>
#include <fstream>
#include <iostream>
#include <cassert>
#include "BigTorus.h"
#include "TRGSW.h"

using namespace std;

void read_lwe_key(const char* const fn, int64_t* const key, const int64_t n) {
    const int32_t LWE_KEY_TYPE_UID = 43;

    ifstream f(fn, ifstream::binary);
    if (not f.is_open()) {
        fprintf(stderr, "Function %s: cannot open file %s\n", __FUNCTION__, fn);
        exit(-1);
    }

    int32_t type_uid;
    f.read((char*)&type_uid, sizeof(int32_t));
    assert(type_uid == LWE_KEY_TYPE_UID);

    int* key_tmp = new int[n];
    f.read((char*)key_tmp, sizeof(int)*n);

    for (int i = 0; i < n; ++i)
        key[i] = (int64_t)key_tmp[i];
    delete[] key_tmp;

    f.close();

    // printf("Read key:\n");
    // for (int i = 0; i < n; ++i)
    //     printf("%ld ", key[i]);
    // printf("\n");
}

int main() {

    // Generation of Bootstrapping key (s)

    int64_t N = 4096;
    int64_t n_in = 500;
    int64_t nblimbs = 2;
    int64_t alpha_bits = 120; //signed

    BigTorusParams bt_params(nblimbs);

    TRGSWParams trgswParams(N, bt_params);

    BigTorusPolynomial phase(N, bt_params);

    TRLwe reps(trgswParams);
    TRGSW *c = new_TRGSW_array(n_in, trgswParams);
    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);

    int64_t *s = new int64_t[n_in];

    cout << "start encrypt bootstrapping key: " << clock() / double(CLOCKS_PER_SEC) << endl;

    read_lwe_key("secret_keyset.bin", s, n_in);

#pragma omp parallel for
    for (int i = 0; i < n_in; i++) {

        int_encrypt(c[i], s[i], *key, alpha_bits);

    }

    cout << "end encrypt bootstrapping key: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //Generation KS key (ks_key)

    int64_t N_in = N;
    int64_t N_out = N;
    int64_t nblimbs_in = 3;
    int64_t nblimbs_out = 3;
    int64_t alpha_bits_ks = 80;
    //int64_t limb_prec = limb_precision(alpha_bits_ks);

    BigTorusParams bt_params_in(nblimbs_in);
    BigTorusParams bt_params_out(nblimbs_out);
    TLweParams tlwe_params_in(N_in, bt_params_in);
    TRLweParams trlwe_params_out(N_out, bt_params_out);

    shared_ptr<TLweKey> key_in = key;
    shared_ptr<TLweKey> key_out = key;

    cout << "start keygen key switch key at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    shared_ptr<pubKsKey32> ks_key = ks_keygen32(
            trlwe_params_out, tlwe_params_in,
            *key_in, *key_out, alpha_bits_ks);
    cout << "end keygen key switch key at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //Generation RK key (rk)

    int64_t alpha_rk = 110;
    int64_t nblimbs_rk = limb_precision(alpha_rk);
    BigTorusParams bt_params_rk(nblimbs_rk);
    TRGSWParams trgswParams_rk(N, bt_params_rk);
    TRGSW rk(trgswParams_rk);
    cout << "start keygen rk key at: " << clock() / double(CLOCKS_PER_SEC) << endl;
    intPoly_encrypt(rk, key->key, *key, alpha_rk);
    cout << "end keygen rk key at: " << clock() / double(CLOCKS_PER_SEC) << endl;

    //serialize bootstrapping key
    ofstream bk_key_out("bk.key");
    ostream_write_binary(bk_key_out, &n_in, sizeof(int64_t));
    serializeTRGSWParams(bk_key_out, trgswParams);
    for (int j = 0; j < n_in; ++j) {
        serializeTRGSWContent(bk_key_out, c[j]);
    }
    bk_key_out.close();

    //serialize key switch key
    ofstream ks_key_out("ks.key");
    serializepubKsKey32(ks_key_out, *ks_key);
    ks_key_out.close();

    //serialize relian. key
    ofstream rk_key_out("rk.key");
    serializeTRGSW(rk_key_out, rk);
    rk_key_out.close();


}
