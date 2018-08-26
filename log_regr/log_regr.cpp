#include "lr_params.h"
#include "keyset.h"
#include "log_regr_fncs.h"
#include "io_ctxt.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "cxxopts.hpp"


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
    compute_X_beta(X_beta_cp, beta_in_out, X_cols, keyset->ks_l2_l1, params->trgsw_params_l1, lr_params);

    // #ifdef DEBUG
    // {
    //     printf("######### DEBUG MSG BEG #########\n");
    //     for (int i = 0; i < lr_params.k; ++i)
    //     {
    //         printf("X_{:,%d}.beta (1xTRLWE): ", i);
    //         // double scale = 1.0;
    //         double scale = 1/lr_params.X_beta_scale;
    //         print_trlwe_sample(X_beta_cp+i, secret_keyset->trlwe_key_l1, 10, scale);
    //         printf("\n");
    //     }
    //     printf("######### DEBUG MSG END #########\n");
    // }
    // #endif


    VERBOSE_1_PRINT("==> accumulate columns of X.beta ... ");
    for (int j = 1; j < lr_params.k; ++j)
        TLweFunctions<Torus>::AddTo(X_beta_cp+0, X_beta_cp+j, params->trlwe_params_l1);
    VERBOSE_1_PRINT("done\n");

    #ifdef DEBUG
    {
        printf("######### DEBUG MSG BEG #########\n");
        printf("X.beta (1xTRLWE): ");
        // double scale = 1.0;
        double scale = 1.0/lr_params.X_beta_scale;
        print_trlwe_sample(X_beta_cp+0, secret_keyset->trlwe_key_l1, 10, scale);
        printf("\n");
        printf("######### DEBUG MSG END #########\n");
    }
    #endif

    LweSample<Torus>* X_beta = nullptr;
    if (X_beta_out != nullptr)
        X_beta = X_beta_out;
    else
        X_beta = new_obj_array<LweSample<Torus>>(lr_params.n, params->tlwe_params_l0);

    extract_X_beta(X_beta, X_beta_cp, params->trlwe_params_l1, keyset->ks_l1_l0, lr_params);
    del_obj_array(lr_params.k, X_beta_cp);

    #ifdef DEBUG
        printf("######### DEBUG MSG BEG #########\n");
        printf("X.beta (%dxTLWE): ", lr_params.n);
        for (int i = 0; i < 10; ++i) {
            print_tlwe_sample(X_beta+i, secret_keyset->tlwe_key_l0, 1./lr_params.X_beta_scale);
        }
        printf("\n");
        printf("######### DEBUG MSG END #########\n");
    #endif

    if (X_beta_out != nullptr) return;

    TLweSample<Torus> *acc = new_obj_array<TLweSample<Torus>>(lr_params.n, params->trlwe_params_l2);
    compute_Xt_sigma(acc, sigmoid_xt_tps, params->trgsw_params_l2, X_beta, params->tlwe_params_l0, keyset->bk_fft, lr_params);

    if (X_beta_out == nullptr)
        del_obj_array(lr_params.n, X_beta);

    VERBOSE_1_PRINT("==> accumulate to compute \\alpha.Xt.\\sigma(X.beta) ... ");
    for (int i = 1; i < lr_params.n; ++i) {
        TLweFunctions<Torus>::AddTo(acc+0, acc+i, params->trlwe_params_l2);
    }
    VERBOSE_1_PRINT("done\n");

    extract_sub_Xt_y(beta_in_out, acc, Xt_y_arr, params->trlwe_params_l2, lr_params);
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

    #ifdef DEBUG
    {
        printf("######### DEBUG MSG BEG #########\n");
        printf("Xt.y: ");
        double scale = 1/lr_params.X_y_scale;
        for (int j = 0; j < lr_params.k; ++j) {
            print_tlwe_sample(Xt_y_arr+j, secret_keyset->tlwe_key_l2, scale);
        }
        printf("\n");
        printf("######### DEBUG MSG END #########\n");
    }
    #endif

    VERBOSE_1_PRINT("Compute initial beta at L2... ");
    LweSample<Torus>* beta = new_obj_array<LweSample<Torus>>(lr_params.k, params->tlwe_params_l2);
    compute_initial_beta(beta, X_cols_l2, params->trgsw_params_l2, y, params->trlwe_params_l2, lr_params);
    VERBOSE_1_PRINT("done\n");

    #ifdef DEBUG
    {
        printf("######### DEBUG MSG BEG #########\n");
        printf("beta^(%d): ", 0);
        // double scale = 1.0;
        double scale = 1/lr_params.beta_scale;
        for (int j = 0; j < lr_params.k; ++j) {
            print_tlwe_sample(beta+j, secret_keyset->tlwe_key_l2, scale);
        }
        printf("\n");
        printf("######### DEBUG MSG END #########\n");
    }
    #endif

    write_tlwe_samples(lr_params.filename_prefix_beta + to_string(0) + ".ctxt", beta, params->tlwe_params_l2, lr_params.k);

    VERBOSE_1_PRINT("Logistic regression:\n");
    for (int iter = 1; iter < lr_params.nb_iters; ++iter)
    {
        VERBOSE_1_PRINT("> iteration %d start\n", iter);
        log_regr_iter(beta, sigmoid_xt_tps, Xt_y_arr, X_cols_l1, keyset, params, lr_params);

        #ifdef DEBUG
            printf("######### DEBUG MSG BEG #########\n");
            printf("beta^(%d): ", iter);
            for (int j = 0; j < lr_params.k; ++j) {
                print_tlwe_sample(beta+j, secret_keyset->tlwe_key_l2, 1/lr_params.beta_scale);
            }
            printf("\n");
            printf("######### DEBUG MSG END #########\n");
        #endif

        write_tlwe_samples(lr_params.filename_prefix_beta + to_string(iter) + ".ctxt", beta, params->tlwe_params_l2, lr_params.k);
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
        // ("inp_file",          "encrypted data input file", cxxopts::value<string>()->default_value(lr_params.filename_data))
        // ("cloud_keyset",      "cloud keyset file", cxxopts::value<string>()->default_value(lr_params.filename_cloud_keyset))
        // ("he_params",         "HE parameters files", cxxopts::value<string>()->default_value(lr_params.filename_params))
        // ("out_path_prefix",   "output path prefix", cxxopts::value<string>()->default_value(""))
        ("iters",             "number of iterations", cxxopts::value<uint>()->default_value("5"))
        ("seed",              "random generator seed", cxxopts::value<uint>()->default_value("42"))
        ;

    options.add_options("general")
        ("threads",           "number of execution threads", cxxopts::value<uint>()->default_value("4"))
        ("v",                 "verbose level ('v' count gives verbosity level)")
        ("h,help",            "print this message")
        ;

    auto result = options.parse(argc, argv);

    if (result.count("h") > 0) {
        printf("%s\n", options.help(options.groups()).c_str());
        exit(-1);
    }

    lr_params.nb_iters = result["iters"].as<uint>();
    lr_params.seed = result["seed"].as<uint>();
    lr_params.nb_threads = result["threads"].as<uint>();
    lr_params.verbose_level = result.count("v");

    return lr_params;
}

int main(int argc, char **argv) {
    LRParams lr_params = parse_args(argc, argv);

    const TfheParamSet *params = TfheParamSet::read(lr_params.filename_params);
    const TfheCloudKeySet *keyset = TfheCloudKeySet::read(lr_params.filename_cloud_keyset, params);

    #ifdef DEBUG
        secret_keyset = TfheSecretKeySet::read(lr_params.filename_secret_keyset, params);
    #endif

    TGswSample<Torus>* X_cols_l1 = nullptr;
    TGswSample<Torus>* X_cols_l2 = nullptr;
    TLweSample<Torus>* y = nullptr;
    TLweSample<Torus>* sigmoid_xt_tps = nullptr;

    read_data(lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);

    if (lr_params.verbose_level)
        lr_params.print();

    RandomGen::set_seed(lr_params.seed);

    #ifdef DEBUG
        printf("######### DEBUG MSG BEG #########\n");

        printf("Input data:\n");

        printf("X_cols_l1:\n");
        for (int j = 0; j < lr_params.k; ++j) {
            printf("col %1d:\n", j);
            print_X_col(X_cols_l1+j, secret_keyset->trgsw_key_l1, lr_params);
            printf("\n");
        }

        // printf("X_cols_l2:\n");
        // for (int j = 0; j < lr_params.k; ++j) {
        //     printf("col %1d:\n", j);
        //     print_trgsw_sample(X_cols_l2+j, secret_keyset->trgsw_key_l2, lr_params.n, lr_params.X_range);
        //     printf("\n");
        // }

        printf("y[::-1]:\n");
        // print_y(y, secret_keyset->trlwe_key_l2, lr_params);
        print_trlwe_sample(y, secret_keyset->trlwe_key_l2, lr_params.n, 1/lr_params.y_scale);
        printf("\n");

        // printf("sigmoid_xt_tps:\n");
        // for (int i = 0; i < lr_params.n; ++i) {
        //     printf("tp %3d:\n", i);
        //     print_trlwe_sample(sigmoid_xt_tps+i, secret_keyset->trlwe_key_l2, params->trlwe_params_l2->N, 1/lr_params.y_scale);
        //     printf("\n");
        // }
        printf("######### DEBUG MSG END #########\n");
    #endif

    log_regr(X_cols_l1, X_cols_l2, y, sigmoid_xt_tps, keyset, params, lr_params);

    del_obj_array(lr_params.n, sigmoid_xt_tps);
    del_obj(y);
    del_obj_array(lr_params.k, X_cols_l2);
    del_obj_array(lr_params.k, X_cols_l1);

    #ifdef DEBUG
        delete secret_keyset;
    #endif

    delete keyset;
    delete params;

    return 0;
}
