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

    int X_precbit = 8;

    //
    int nb_iters;
    int nb_threads;
    int seed;
    int verbose_level;

    int sigm_levels = 16;
    double sigm_in_min = -2.0;
    double sigm_in_max = 2.0;
    double scale = 1.0 / 2;

    const char* const params_filename = "params.bin";
    const char* const secret_keyset_filename = "secret_keyset.bin";
    const char* const cloud_keyset_filename = "cloud_keyset.bin";

    const char* const filename_data = "data.ctxt";

    const char* const filename_X_beta = "X_beta.ctxt";
    const char* const filename_prefix_beta = "beta_iter_";

    /* distance between k consecutive coefficients encoded in test polynomials */
    int d;

    /* LR step */
    const int alpha = 4;

    /**
     * @brief Update computed options
     */
    void update(int N = 8192) {
        d = (int)((float)N  / k / sigm_levels);
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
