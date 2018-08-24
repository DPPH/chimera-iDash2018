#ifndef IO_CTXT_H
#define IO_CTXT_H

#include "tfhe_core.h"
#include "tfhe_io.h"

#include <iostream>
#include <fstream>

using namespace std;

// void write_trgsw_samples(const string& filename, const TGswSample<Torus>* samples, const TGswParams<Torus>* params, const int32_t n) {
    // ofstream out(filename, ofstream::binary);
    // StdOstream out_stream = to_Ostream(out);
    // out_stream.fwrite(&n, sizeof(int32_t));
    // for (int i = 0; i < n; ++i)
        // IOFunctions<Torus>::write_lweSample(out_stream, samples+i, params);
    // out.close();
// }

void write_tlwe_samples(const string& filename, const LweSample<Torus>* samples, const LweParams<Torus>* params, const int32_t n) {
    ofstream out(filename, ofstream::binary);
    StdOstream out_stream = to_Ostream(out);
    for (int i = 0; i < n; ++i)
        IOFunctions<Torus>::write_lweSample(out_stream, samples+i, params);
    out.close();
}

LweSample<Torus>* read_tlwe_samples(const string& filename, const LweParams<Torus>* params, const int32_t n) {
    ifstream inp(filename, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    LweSample<Torus>* samples = new_obj_array<LweSample<Torus>>(n, params);
    for (int i = 0; i < n; ++i)
        IOFunctions<Torus>::read_lweSample(inp_stream, samples+i, params);
    inp.close();
    return samples;
}


TGswSample<Torus>* read_trgsw_samples(StdIstream& inp_stream, const TGswParams<Torus> *params, const int32_t n) {
    TGswSample<Torus>* samples = new_obj_array<TGswSample<Torus>>(n, params);
    for (int j = 0; j < n; ++j)
        IOFunctions<Torus>::read_tGswSample(inp_stream, samples+j, params);
    return samples;
}

TLweSample<Torus>* read_trlwe_samples(StdIstream& inp_stream, const TLweParams<Torus> *params, const int32_t n) {
    TLweSample<Torus> *samples = new_obj_array<TLweSample<Torus>>(n, params);
    for (int i = 0; i < n; ++i)
        IOFunctions<Torus>::read_tLweSample(inp_stream, samples+i, params);
    return samples;
}

void read_data_header(
    StdIstream& inp_stream,
    LRParams& lr_params
    )
{
    inp_stream.fread(&lr_params.n, sizeof(int32_t));
    inp_stream.fread(&lr_params.m, sizeof(int32_t));
    inp_stream.fread(&lr_params.k, sizeof(int32_t));

    char buf[18];
    inp_stream.fread(buf, 18);
    sscanf(buf, "%f", &(lr_params.X_scale));
    // printf("X scale %s %f\n", buf, lr_params.X_scale);

    lr_params.update();
}

void read_data_header(
    LRParams& lr_params
    )
{
    ifstream inp(lr_params.filename_data, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    read_data_header(inp_stream, lr_params);
    inp.close();
}

void read_data(
    StdIstream& inp_stream,
    LRParams& lr_params,
    TLweSample<Torus>*& sigmoid_xt_tps,
    TLweSample<Torus>*& y,
    TGswSample<Torus>*& X_cols_l1,
    TGswSample<Torus>*& X_cols_l2,
    const TfheParamSet* params)
{
    read_data_header(inp_stream, lr_params);

    /* read sigmoid . Xt test polynomials */
    sigmoid_xt_tps = read_trlwe_samples(inp_stream, params->trlwe_params_l2, lr_params.n);

    /* read y */
    y = read_trlwe_samples(inp_stream, params->trlwe_params_l2, 1);

    /* read X_cols L1 */
    X_cols_l1 = read_trgsw_samples(inp_stream, params->trgsw_params_l1, lr_params.k);

    /* read X_cols L2 */
    X_cols_l2 = read_trgsw_samples(inp_stream, params->trgsw_params_l2, lr_params.k);
}

void read_data(
    LRParams& lr_params,
    TLweSample<Torus>*& sigmoid_xt_tps,
    TLweSample<Torus>*& y,
    TGswSample<Torus>*& X_cols_l1,
    TGswSample<Torus>*& X_cols_l2,
    const TfheParamSet* params)
{
    ifstream inp(lr_params.filename_data, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    read_data(inp_stream, lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);
    inp.close();
}


void print_trgsw_sample(const TGswSample<Torus>* sample, const TGswKey<Torus>* key, const int nb_coefs, const int Msize) {
    const TGswParams<Torus>* params = key->params;
    IntPolynomial* result = new_obj<IntPolynomial>(params->tlwe_params->N);
    TGswFunctions<Torus>::SymDecrypt(result, sample, key, Msize);

    for (int i = 0; i < nb_coefs; ++i) {
        if (result->coefs[i] <= Msize/2)
            printf("%3d ", result->coefs[i]);
        else
            printf("%3d ", result->coefs[i]-Msize);
    }
    del_obj(result);
}

void print_trlwe_sample(const TLweSample<Torus>* sample, const TLweKey<Torus>* key, const int nb_coefs, const float scale=1.0) {
    TorusPolynomial<Torus>* result = nullptr;

    if (key != nullptr) {
        const TLweParams<Torus>* params = key->params;
        result = new_obj<TorusPolynomial<Torus>>(params->N);
        TLweFunctions<Torus>::Phase(result, sample, key);
    } else {
        result = sample->b;
    }

    for (int i = 0; i < nb_coefs; ++i) {
        printf("%.7lf ", TorusUtils<Torus>::to_double(result->coefsT[i]) * scale);
    }

    if (key != nullptr)
        del_obj(result);
}

void print_trlwe_sample_coef(const TLweSample<Torus>* sample, const TLweKey<Torus>* key, const int idx, const float scale=1.0) {
    TorusPolynomial<Torus>* result = nullptr;

    if (key != nullptr) {
        const TLweParams<Torus>* params = key->params;
        result = new_obj<TorusPolynomial<Torus>>(params->N);
        TLweFunctions<Torus>::Phase(result, sample, key);
    } else {
        result = sample->b;
    }

    printf("%.12lf ", TorusUtils<Torus>::to_double(result->coefsT[idx]) * scale);

    if (key != nullptr)
        del_obj(result);
}

void print_tlwe_sample(const LweSample<Torus>* sample, const LweKey<Torus>* key, const double scale = 1.0) {
    Torus result = LweFunctions<Torus>::Phase(sample, key);
    printf("%.7lf ", TorusUtils<Torus>::to_double(result) * scale);
}


void print_X_col(const TGswSample<Torus>* X_col, const TGswKey<Torus>* key, const LRParams& lr_params, int nb_coefs = -1) {
    if (nb_coefs == -1)
        nb_coefs = lr_params.n;
    int Msize = lr_params.X_range;

    const TGswParams<Torus>* params = key->params;
    IntPolynomial* result = new_obj<IntPolynomial>(params->tlwe_params->N);
    TGswFunctions<Torus>::SymDecrypt(result, X_col, key, Msize);

    for (int i = 0; i < nb_coefs; ++i) {
        int t = result->coefs[i];
        if (result->coefs[i] > Msize/2)
            t -= Msize;
        printf("%.7f ", t / lr_params.X_scale);

    }
    del_obj(result);
}

#endif
