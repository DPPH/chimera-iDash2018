#include <iostream>
#include <fstream>

#include "tfhe_core.h"
#include "tfhe_io.h"
#include "tfhe_garbage_collector.h"

using namespace std;

constexpr inline double mulBySqrtTwoOverPi(double x) { return x*sqrt(2./M_PI); }

TFheGateBootstrappingParameterSet<Torus>* new_bootstrapping_parameters() {
    // const int N = 1024;
    const int N = 1024*16;
    const int k = 1;
    const int n = 500;

    const int bk_l = 8;
    const int bk_Bgbit = 6;

    const int ks_t = 24;
    const int ks_basebit = 1;

    const double ks_stdev = pow(2.,-33);   //standard deviation
    const double bk_stdev = pow(2.,-50);          //standard deviation
    const double max_stdev = mulBySqrtTwoOverPi(pow(2.,-8)/4.);

    LweParams<Torus>* params_in_out = new_obj<LweParams<Torus>>(n, ks_stdev, max_stdev);
    TLweParams<Torus>* params_accum = new_obj<TLweParams<Torus>>(N, k, bk_stdev, max_stdev);
    TGswParams<Torus>* params_bk = new_obj<TGswParams<Torus>>(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector<Torus>::register_param(params_in_out);
    TfheGarbageCollector<Torus>::register_param(params_accum);
    TfheGarbageCollector<Torus>::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet<Torus>(ks_t, ks_basebit, params_in_out, params_bk);
}

int main(int argc, char const *argv[]) {
  uint32_t values[2];
  values[0] = time(NULL)>>32;
  values[1] = time(NULL);
  RandomGen::set_seed(values, 2);

  // generate params
  TFheGateBootstrappingParameterSet<Torus>* params = new_bootstrapping_parameters();

  // generate the secret keyset and obtain cloud keyset
  TFheGateBootstrappingSecretKeySet<Torus>* secret_keyset = TFheGateBootstrappingSecretKeySet<Torus>::new_random(params);
  const TFheGateBootstrappingCloudKeySet<Torus>* cloud_keyset = &(secret_keyset->cloud);

  ofstream f1("secret.key", ofstream::binary);
  export_tfheGateBootstrappingSecretKeySet_toStream(f1, secret_keyset, true);
  ofstream f2("cloud.key", ofstream::binary);
  IOFunctions<Torus>::write_tfheGateBootstrappingCloudKeySet(to_Ostream(f2), cloud_keyset, true);
}
