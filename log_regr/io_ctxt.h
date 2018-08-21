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

#endif
