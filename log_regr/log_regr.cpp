#include "lr_params.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "cxxopts.hpp"

/**
 * @brief Blind rotate \c acc by TLWE sample \c input coefficients
 */
void blindrotate(
    TLweSample<Torus>* acc,
    const LweSample<Torus>* input,
    const LweBootstrappingKeyFFT<Torus>* bk
    )
{
    const TLweParams<Torus>* accum_params = bk->accum_params;
    const int N=accum_params->N;
    const int _2N = 2*N;
    const int n = bk->in_out_params->n;

    // Modulus switching
    int* bara = new int[N];
    int barb = TorusUtils<Torus>::modSwitchFromTorus(input->b,_2N);
    for (int i=0; i<n; i++) {
        bara[i]=TorusUtils<Torus>::modSwitchFromTorus(input->a[i],_2N);
    }

    TLweSample<Torus>* tmp = new_obj<TLweSample<Torus>>(accum_params);
    tLweCopy(tmp, acc, accum_params);

    // Accumulator rotate by -barb
    TLweFunctions<Torus>::MulByXai(acc, _2N-barb, tmp, accum_params);

    // Accumulator blind rotate by \sum bara_i * sk_i
    tfhe_blindRotate_FFT(acc, bk->bkFFT, bara, n, bk->bk_params);

    delete[] bara;
}


void compute_Xt_sigma(
    TLweSample<Torus> *result,
    const LweSample<Torus>* inputs,
    const TLweSample<Torus> *sigmoid_xt_tps,
    const int tps_cnt,
    const LweBootstrappingKeyFFT<Torus>* bk
    )
{
    const TLweParams<Torus> *accum_params = bk->accum_params;

    TLweSample<Torus>* acc = new_obj<TLweSample<Torus>>(accum_params);
    for (int i = 0; i < tps_cnt; ++i) {
        tLweCopy(acc, sigmoid_xt_tps+i, accum_params);
        blindrotate(acc, inputs+i, bk);

        tLweAddTo(result, acc, accum_params);
    }
    del_obj(acc);
}

// void compute_gradient(
//     const TLweSample<Torus> *Xt_y,
//     const LweBootstrappingKeyFFT<Torus>* bk
//     )
// {
//     const TLweParams<Torus> *trlwe_params = bk->accum_params->tlwe_params;
//     TLweSample<Torus>* result = new_obj<TLweSample<Torus>>(trlwe_params);
//     tLweCopy(result, Xt_y, trlwe_params);



//     compute_Xt_sigma()
// }


/**
 * @brief Extract coefficients i*delta for i=0..nb_coefs-1 from \c input TRLWE sample.
 */
void extract_coefs(
    LweSample<Torus> *outputs,
    const TLweSample<Torus> *input,
    const int nb_coefs,
    const int delta
    )
{
    // void TLweFunctions<TORUS>::ExtractLweSampleIndex(LweSample<TORUS>* result, const TLweSample<TORUS>* x,
    //     const int index, const LweParams<TORUS>* params,  const TLweParams<TORUS>* rparams) {

    // void lweKeySwitch(LweSample<TORUS>* result, const LweKeySwitchKey<TORUS>* ks, const LweSample<TORUS>* sample) {


    for (int i = 0; i < nb_coefs; ++i) {
        TLweFunctions<TORUS>::ExtractLweSampleIndex(outputs+i, input, i*delta, ., .);


    }
}


LRParams parse_args(int argc, char** argv) {
    LRParams lr_params;

    cxxopts::Options options("LR", "Logistic regression learning on encrypted genomic data");

    options.add_options("input data")
        ("n,indivs",          "number of individuals", cxxopts::value<uint>()->default_value("245"))
        ("k,feats",           "number of features", cxxopts::value<uint>()->default_value("3"))
        ("m,snps",            "number of SNPs", cxxopts::value<uint>()->default_value("10643"))
        // ("b,batch",           "batch size", cxxopts::value<uint>()->default_value("2"))
        ("iters",             "number of iterations", cxxopts::value<uint>()->default_value("5"))
        ("seed",              "random generator seed", cxxopts::value<uint>()->default_value("42"))
        ;

    options.add_options("general")
        // ("path_train",        "encrypted train data input path", cxxopts::value<string>()->default_value("/dev/shm/enc_data"))
        // ("file_train",        "train input data .csv file (perform learning on clear data)", cxxopts::value<string>())
        // ("file_valid",        "validation input data .csv file (secret key file 'secret.key' must be present)", cxxopts::value<string>())
        // ("path_beta_out",     "beta output path", cxxopts::value<string>()->default_value(""))
        ("threads",           "number of execution threads", cxxopts::value<uint>()->default_value("4"))
        ("v",                 "verbose level ('v' count gives verbosity level)")
        ("h,help",            "print this message")
        ;

    auto result = options.parse(argc, argv);

    if (result.count("h") > 0) {
        printf("%s\n", options.help(options.groups()).c_str());
        exit(-1);
    }

  // if (result.count("file_train") > 0 and result["file_train"].as<string>().length()) {
  //   PARAMS.file_train = result["file_train"].as<string>();
  //   PARAMS.learn_on_clear = true;
  // }
  // else if (result.count("path_train") > 0 and result["path_train"].as<string>().length()) {
  //   PARAMS.path_train = result["path_train"].as<string>();
  //   PARAMS.learn_on_clear = false;
  // } else {
  //   fprintf(stderr, "Please specify input path!!!\n");
  //   fprintf(stderr, "%s\n", result.help(result.groups()).c_str());
  //   exit(-1);
  // }

  // if (result.count("file_valid") > 0 and result["file_valid"].as<string>().length()) {
  //   PARAMS.file_valid = result["file_valid"].as<string>();
  //   PARAMS.valid = true;
  // }


  // if (not PARAMS.learn_on_clear) {
  //   if (result.count("path_beta_out") > 0 and result["path_beta_out"].as<string>().length()) {
  //     PARAMS.path_beta_out = result["path_beta_out"].as<string>();
  //   } else {
  //     fprintf(stderr, "Please specify output theta path!!!\n");
  //     fprintf(stderr, "%s\n", result.help(result.groups()).c_str());
  //     exit(-1);
  //   }
  // }

    lr_params.n = result["indivs"].as<uint>();
    lr_params.k = result["feats"].as<uint>()+1;
    lr_params.m = result["snps"].as<uint>();

    lr_params.nb_iters = result["iters"].as<uint>();
    lr_params.seed = result["seed"].as<uint>();
    lr_params.nb_threads = result["threads"].as<uint>();
    lr_params.verbose_level = result.count("v");

    if (lr_params.verbose_level)
        lr_params.print();

    return lr_params;
}


int main(int argc, char **argv) {
  LRParams lr_params = parse_args(argc, argv);
}
