#include "lr_params.h"
#include "keyset.h"
#include "log_regr_fncs.h"
#include "io_ctxt.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "cxxopts.hpp"

#define DEBUG

void log_regr_iter(
    LweSample<Torus> *beta_in_out, // k TLWE L2 encryptions of \beta
    const TLweSample<Torus> *sigmoid_xt_tps, // TRLWE L2 samples encrypting test polynomials Xt . \sigma(u) for u=-A..A
    const LweSample<Torus>* Xt_y_arr, // k TLWE L2 encryptions of vector Xt.y.\alpha elements
    const TGswSample<Torus>* X_cols, // k TRGSW L1 samples of X columns \sum_{i=0}^{n-1} X_ij Z^i for j=0..k-1

    const TfheCloudKeySet* keyset,
    const TfheParamSet* params,
    const LRParams& lr_params,
    LweSample<Torus>* X_beta_out = nullptr
    )
{
    TLweSample<Torus>* X_beta_cp = new_obj_array<TLweSample<Torus>>(lr_params.k, params->trlwe_params_l1);

    compute_X_beta(X_beta_cp, beta_in_out, X_cols, params->trgsw_params_l1, lr_params);

    VERBOSE_1_PRINT("==> accumulate columns of X.beta ... ");
    for (int j = 1; j < lr_params.k; ++j)
        TLweFunctions<Torus>::AddTo(X_beta_cp+0, X_beta_cp+j, params->trlwe_params_l1);
    VERBOSE_1_PRINT("done\n");

    LweSample<Torus>* X_beta = nullptr;
    if (X_beta_out != nullptr)
        X_beta = X_beta_out;
    else
        X_beta = new_obj_array<LweSample<Torus>>(lr_params.n, params->tlwe_params_l0);

    extract_X_beta(X_beta, X_beta_cp, params->trlwe_params_l1, keyset->ks_l1_l0, lr_params);
    del_obj_array(lr_params.k, X_beta_cp);

    if (X_beta_out != nullptr) return;

    TLweSample<Torus> *acc = new_obj_array<TLweSample<Torus>>(lr_params.n, params->trlwe_params_l2);
    compute_Xt_sigma(acc, sigmoid_xt_tps, params->trgsw_params_l2, X_beta, params->tlwe_params_l0, keyset->bk_fft, lr_params);

    if (X_beta_out == nullptr)
        del_obj_array(lr_params.n, X_beta);

    accumulate_trlwe_acc(acc, params->trlwe_params_l2, lr_params);

    extract_add_Xt_y(beta_in_out, acc, Xt_y_arr, params->trlwe_params_l2, lr_params);
    del_obj_array(lr_params.n, acc);
}

void log_regr(
    const TGswSample<Torus>* X_cols_l1,
    const TGswSample<Torus>* X_cols_l2,
    const TLweSample<Torus>* y,
    const TLweSample<Torus> *sigmoid_xt_tps,

    const TfheCloudKeySet* keyset,
    const TfheParamSet* params,
    const LRParams& lr_params
    )
{
    /* pre-compute Xt.y */
    VERBOSE_1_PRINT("Compute Xt.y at L2... ");
    LweSample<Torus>* Xt_y_arr = new_obj_array<LweSample<Torus>>(lr_params.k, params->tlwe_params_l2);
    compute_Xt_y(Xt_y_arr, X_cols_l2, y, params->trgsw_params_l2, lr_params);
    VERBOSE_1_PRINT("done\n");

    VERBOSE_1_PRINT("Compute initial beta at L2... ");
    LweSample<Torus>* beta = new_obj_array<LweSample<Torus>>(lr_params.k, params->tlwe_params_l2);
    compute_initial_beta(beta, X_cols_l2, params->trgsw_params_l2, y, params->trlwe_params_l2, lr_params);
    VERBOSE_1_PRINT("done\n");

    write_tlwe_samples(lr_params.filename_prefix_beta + to_string(0) + ".ctxt", beta, params->tlwe_params_l1, lr_params.k);

    VERBOSE_1_PRINT("Logistic regression:\n");
    for (int iter = 1; iter < lr_params.nb_iters; ++iter)
    {
        VERBOSE_1_PRINT("> iteration %d start\n", iter);
        log_regr_iter(beta, sigmoid_xt_tps, Xt_y_arr, X_cols_l1, keyset, params, lr_params);
        write_tlwe_samples(lr_params.filename_prefix_beta + to_string(iter) + ".ctxt", beta, params->tlwe_params_l1, lr_params.k);
    }

    VERBOSE_1_PRINT("> compute X.beta result\n");
    LweSample<Torus>* X_beta = new_obj_array<LweSample<Torus>>(lr_params.n, params->tlwe_params_l0);
    log_regr_iter(beta, sigmoid_xt_tps, Xt_y_arr, X_cols_l1, keyset, params, lr_params, X_beta);

    VERBOSE_1_PRINT("> write X.beta\n");
    write_tlwe_samples(lr_params.filename_X_beta, X_beta, params->tlwe_params_l0, lr_params.n);

    del_obj_array(lr_params.n, X_beta);
    del_obj_array(lr_params.k, Xt_y_arr);
    del_obj_array(lr_params.k, beta);
}

LRParams parse_args(int argc, char** argv) {
    LRParams lr_params;

    cxxopts::Options options("LR", "Logistic regression learning on encrypted genomic data");

    options.add_options("input data")
        ("n,indivs",          "number of individuals", cxxopts::value<uint>()->default_value("245"))
        ("k,feats",           "number of features", cxxopts::value<uint>()->default_value("3"))
        ("m,snps",            "number of SNPs", cxxopts::value<uint>()->default_value("10643"))
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

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheCloudKeySet *keyset = TfheCloudKeySet::read(lr_params.cloud_keyset_filename, params);

    TGswSample<Torus>* X_cols_l1 = nullptr;
    TGswSample<Torus>* X_cols_l2 = nullptr;
    TLweSample<Torus>* y = nullptr;
    TLweSample<Torus>* sigmoid_xt_tps = nullptr;

    read_data(lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);

    lr_params.n = 16;
    lr_params.nb_iters = 2;
    lr_params.update();

    log_regr(X_cols_l1, X_cols_l2, y, sigmoid_xt_tps, keyset, params, lr_params);

    return 0;
}
