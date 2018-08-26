#ifndef LRPARAMS_H
#define LRPARAMS_H

#include <cmath>
#include <stdio.h>

struct LRParams
{
    // data set shape
    int n;          // number X rows
    int k;          // number X cols
    int m;          // number SNPs

    const int X_range = 1<<8; // 8 bits for X matrix
    float X_scale;
    float y_scale;
    float sigmoid_scale;
    float X_beta_scale;
    float X_y_scale;
    float beta_scale;

    // const int y_precbit = 19;
    // const double y_scale = pow(2., -y_precbit);
    // const int y_range = 1<<y_precbit;

    //
    int nb_iters;
    int nb_threads;
    int seed;
    int verbose_level;

    // const int sigm_levels = 16;
    const double sigm_in_min = -4.0;
    const double sigm_in_max = 4.0;
    // const double sigm_multiplier = pow(2., -6);

    const char* const params_filename = "params.bin";
    const char* const secret_keyset_filename = "secret_keyset.bin";
    const char* const cloud_keyset_filename = "cloud_keyset.bin";
    const char* const ks_l2_l1_filename = "ks_l2_l1.bin";

    const char* const filename_data = "data.ctxt";

    const char* const filename_X_beta = "X_beta.ctxt";
    const char* const filename_prefix_beta = "beta_iter_";

    /* LR step */
    const int alpha = 4;

    /**
     * @brief Update computed options
     */
    void update() {
        y_scale = float(alpha) / 16 / X_scale / X_scale;
        X_beta_scale = 1. / 16;
        sigmoid_scale = X_y_scale = y_scale * X_scale;
        beta_scale = 1. / 16 / X_scale;
    }

    void print() {
        printf("LRParams\n");

        printf("  data:\n");
        printf("    %d individuals\n", n);
        printf("    %d features\n", k-1);
        printf("    %d SNPs\n", m);

        printf("  algorithm:\n");
        printf("    %d iterations\n", nb_iters);
        printf("    %d threads\n", nb_threads);
        printf("    %d seed\n", seed);
        printf("    %d verbose level\n", verbose_level);

        printf("  config:\n");
        printf("    alpha =         %12d\n", alpha);
        printf("    X_scale =       %12.6f\n", X_scale);
        printf("    y_scale =       %12.8f\n", y_scale);
        printf("    X_beta_scale =  %12.8f\n", X_beta_scale);
        printf("    X_y_scale =     %12.8f\n", X_y_scale);
        printf("    beta_scale =    %12.8f\n", beta_scale);
    }
};


#define VERBOSE_1_PRINT(...) \
    if (lr_params.verbose_level > 0) { \
        printf(__VA_ARGS__); \
        fflush(stdout); \
    }
#define VERBOSE_2_PRINT(...) \
    if (lr_params.verbose_level > 1) { \
        printf(__VA_ARGS__); \
        fflush(stdout); \
    }


#endif
