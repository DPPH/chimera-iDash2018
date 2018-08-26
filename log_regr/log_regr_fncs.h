#include "lr_params.h"
#include "keyset.h"
#include "io_ctxt.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "tlwe-functions-extra.h"

#define DEBUG

#ifdef DEBUG
    const TfheSecretKeySet* secret_keyset = nullptr;


    /**
     * @brief Blind rotate \c acc by TLWE sample \c input coefficients
     */
    void blindrotate_fake(
        TLweSample<Torus>* acc,
        const TLweParams<Torus>* acc_params,
        const LweSample<Torus>* inp,
        const LweParams<Torus>* inp_params,
        const TGswSampleFFT<Torus>* bk_fft,
        const TGswParams<Torus>* bk_params,
        const LRParams& lr_params
        )
    {
        const int N=acc_params->N;
        const int _2N = 2*N;
        const int n = inp_params->n;

        // Modulus switching
        int* bara = new int[N];
        int barb = TorusUtils<Torus>::modSwitchFromTorus(inp->b,_2N/4)*4;
        for (int i=0; i<n; i++) {
            bara[i]=TorusUtils<Torus>::modSwitchFromTorus(inp->a[i],_2N/4)*4;
        }

        TorusPolynomial<Torus>* tmp = new_obj<TorusPolynomial<Torus>>(N);
        TLweFunctions<Torus>::Phase(tmp, acc, secret_keyset->trlwe_key_l2);

        // Accumulator rotate by -barb
        TorusPolynomial<Torus>* acc_decr = new_obj<TorusPolynomial<Torus>>(N);
        if (barb != 0)
            TorusPolyFunctions<Torus>::MulByXai(acc_decr, _2N-barb, tmp);

        for (int i = 0; i < n; ++i) {
            if (secret_keyset->tlwe_key_l0->key[i]) {
                TorusPolyFunctions<Torus>::Copy(tmp, acc_decr);
                TorusPolyFunctions<Torus>::MulByXai(acc_decr, bara[i], tmp);
            }
        }
        del_obj(tmp);

        TLweFunctions<Torus>::SymEncrypt(acc, acc_decr, secret_keyset->trlwe_key_l2->params->alpha_min, secret_keyset->trlwe_key_l2);
        del_obj(acc_decr);

        delete[] bara;
    }

    //keyswitch using re-encryption, fast way :)
    void TLweKeySwitch_fake(
        TLweSample<Torus>* out,
        const void* ks_l2_l1,
        const LweSample<Torus>* inp
        )
    {
        const TLweKey<Torus>* out_key = secret_keyset->trlwe_key_l1;
        const LweKey<Torus>* inp_key = secret_keyset->tlwe_key_l2;

        Torus msg_inp = LweFunctions<Torus>::Phase(inp, inp_key);
        TLweFunctions<Torus>::SymEncryptT(out, msg_inp, pow(2,-48), out_key);
    }

    //keyswitch using re-encryption, fast way :)
    void LweKeySwitch_fake(
        LweSample<Torus>* out,
        const void* ks_l1_l0,
        const LweSample<Torus>* inp
        )
    {
        const LweKey<Torus>* out_key = secret_keyset->tlwe_key_l0;
        const LweKey<Torus>* inp_key = secret_keyset->tlwe_key_l1;

        Torus msg_inp = LweFunctions<Torus>::Phase(inp, inp_key);
        LweFunctions<Torus>::SymEncrypt(out, msg_inp, pow(2,-48), out_key);
    }


#endif

/**
 * @brief Blind rotate \c acc by TLWE sample \c input coefficients
 */
void blindrotate(
    TLweSample<Torus>* acc,
    const TLweParams<Torus>* acc_params,
    const LweSample<Torus>* inp,
    const LweParams<Torus>* inp_params,
    const TGswSampleFFT<Torus>* bk_fft,
    const TGswParams<Torus>* bk_params,
    const LRParams& lr_params
    )
{
    const int N=acc_params->N;
    const int _2N = 2*N;
    const int n = inp_params->n;

    // Modulus switching
    int* bara = new int[N];
    int barb = TorusUtils<Torus>::modSwitchFromTorus(inp->b,_2N/4)*4;
    for (int i=0; i<n; i++) {
        bara[i]=TorusUtils<Torus>::modSwitchFromTorus(inp->a[i],_2N/4)*4;
    }

    TLweSample<Torus>* tmp = new_obj<TLweSample<Torus>>(acc_params);
    tLweCopy(tmp, acc, acc_params);

    // Accumulator rotate by -barb
    if (barb != 0)
        TLweFunctions<Torus>::MulByXai(acc, _2N-barb, tmp, acc_params);

    // Accumulator blind rotate by \sum bara_i * sk_i
    tfhe_blindRotate_FFT(acc, bk_fft, bara, n, bk_params);

    delete[] bara;
}

// #define blindrotate blindrotate_fake
#define blindrotate blindrotate

// #define TLweKeySwitch TLweKeySwitch_fake
#define TLweKeySwitch TLweFunctionsExtra<Torus>::KeySwitch

// #define LweKeySwitch LweKeySwitch_fake
#define LweKeySwitch LweFunctions<Torus>::KeySwitch

/**
 * @brief Compute Xt.y.\alpha from TRGSW samples of X columns and TRLWE sample of y (coefficient packed)
 *
 * @param Xt_y k output TLWE samples under extracted TLWE key from TRGSW sample
 * @param X_cols k TRGSW samples encoding \sum_{i=0..n-1} Z^i . X_ij for j=0..k
 * @param y TRLWE sample encodin \sum_{i=0..n-1} Z^(n-i) . y_i
 * @param trgsw_params parameters of TRGSW sample
 * @param n number of rows in matrix X
 * @param k number of columns in matrix X
 */
void compute_Xt_y(
    LweSample<Torus>* Xt_y,
    const TGswSample<Torus>* X_cols,
    const TLweSample<Torus>* y,
    const TGswParams<Torus>* trgsw_params,
    const LRParams& lr_params
    )
{
    const TLweParams<Torus>* trlwe_params = trgsw_params->tlwe_params;
    const LweParams<Torus>* tlwe_params = &(trlwe_params->extracted_lweparams);

    for (int i = 0; i < lr_params.k; ++i) {
        TLweSample<Torus>* acc = new_obj<TLweSample<Torus>>(trlwe_params);
        TGswFunctions<Torus>::ExternProduct(acc, X_cols+i, y, trgsw_params);

        //extract the n-th coefficient
        TLweFunctions<Torus>::ExtractLweSampleIndex(Xt_y+i, acc, lr_params.n-1, tlwe_params, trlwe_params);

        del_obj(acc);
    }
}


/**
 * @brief Compute Xt.(1/2-y)
 */
void compute_initial_beta(
    LweSample<Torus>* beta,
    const TGswSample<Torus>* X_cols,
    const TGswParams<Torus>* trgsw_params,
    const TLweSample<Torus>* y,
    const TLweParams<Torus>* trlwe_params,
    const LRParams& lr_params
    )
{
    /* encrypt 1/2 in first n coefficients */
    TorusPolynomial<Torus>* prob_poly = new_obj<TorusPolynomial<Torus>>(trlwe_params->N);
    TorusPolyFunctions<Torus>::Clear(prob_poly);
    Torus mu = TorusUtils<Torus>::from_double(1./2 * lr_params.y_scale);
    for (int i = 0; i < lr_params.n; ++i)
        prob_poly->coefsT[i] = mu;

    TLweSample<Torus>* ymp = new_obj<TLweSample<Torus>>(trlwe_params);
    TLweFunctions<Torus>::Copy(ymp, y, trlwe_params);
    TorusPolyFunctions<Torus>::SubTo(ymp->b, prob_poly);
    del_obj(prob_poly);

    compute_Xt_y(beta, X_cols, ymp, trgsw_params, lr_params);

    del_obj(ymp);
}


/**
 * @brief Switch \c beta_in TLWE sample to TRLWE and multiply with columns of X
 *
 * @param X_beta_cp TRLWE sample array (k elements), outputs X_{,j} . \beta_j
 * @param beta_in TLWE sample array (k elements) of betas
 * @param X_cols TRGSW sample array (k elements) of X columns
 * @param trgsw_params TRGSW parameters
 */
void compute_X_beta(
    TLweSample<Torus>* X_beta_cp,
    const LweSample<Torus> *beta_in,
    const TGswSample<Torus>* X_cols,
    const TLweKeySwitchKey<Torus>* ks_l2_l1,
    const TGswParams<Torus>* trgsw_params,
    const LRParams& lr_params
    )
{
    /* switch beta from TLWE to TRLWE and multiply each beta by X columns */
    VERBOSE_1_PRINT("==> compute individual columns of X.beta\n");
    #pragma omp parallel for ordered schedule(static,1)
    for (int j = 0; j < lr_params.k; ++j) {
        /* switch from L2 TLWE sample to L1 TRLWE sample  */
        TLweKeySwitch(X_beta_cp+j, ks_l2_l1, beta_in+j);

        // #ifdef DEBUG
        //     printf("######### DEBUG MSG BEG #########\n");
        //     printf("beta_in_%d TLWE: ", j);
        //     print_tlwe_sample(beta_in+j, secret_keyset->tlwe_key_l2, 1./lr_params.beta_scale);
        //     printf("\n");

        //     printf("beta_in_%d TRLWE: ");
        //     print_trlwe_sample(X_beta_cp+j, secret_keyset->trlwe_key_l1, 3, 1./lr_params.beta_scale);
        //     printf("\n");
        //     printf("######### DEBUG MSG END #########\n");
        // #endif


        TGswFunctions<Torus>::ExternMulToTLwe(X_beta_cp+j, X_cols+j, trgsw_params);

        // #ifdef DEBUG
        //     printf("######### DEBUG MSG BEG #########\n");
        //     printf("X_cols: ");
        //     print_trgsw_sample(X_cols+j, secret_keyset->trgsw_key_l1, 10, 1<<lr_params.y_precbit);
        //     printf("\n");

        //     printf("X.beta: ");
        //     print_trlwe_sample(X_beta_cp+j, secret_keyset->trlwe_key_l1, 10);
        //     printf("\n");
        //     printf("######### DEBUG MSG END #########\n");
        // #endif

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", j+1, lr_params.k);
    }
    VERBOSE_2_PRINT("\n");
}



/**
 * @brief Extract individual TLWE samples from TRLWE sample coefficients and keyswitch
 *
 * @param X_beta TLWE sample array (n elements) for output
 * @param X_beta_cp TRLWE sample whose first n coefficients are to be extracted
 */
void extract_X_beta(
    LweSample<Torus>* X_beta,
    const TLweSample<Torus>* X_beta_cp,
    const TLweParams<Torus>* trlwe_params,
    const LweKeySwitchKey<Torus>* ks_l1_l0,
    const LRParams& lr_params
    )
{
    const LweParams<Torus>* tlwe_params = &(trlwe_params->extracted_lweparams);

    /* extract individual TLWE L0 encryptions of \sum_{j=0..k-1} X_ij.\beta_j from X_beta_cp TRLWE L1 sample */
    VERBOSE_1_PRINT("==> extract individual samples from X.beta\n");
    #pragma omp parallel for ordered schedule(static,1)
    for (int i = 0; i < lr_params.n; ++i) {
        LweSample<Torus>* tmp = new_obj<LweSample<Torus>>(tlwe_params);
        TLweFunctions<Torus>::ExtractLweSampleIndex(tmp, X_beta_cp, i, tlwe_params, trlwe_params);

        // #ifdef DEBUG
        // {
        //     // write_tlwe_samples("ks_inp.ctxt", tmp, secret_keyset->tlwe_key_l1->params, 1);

        //     // printf("######### DEBUG MSG BEG #########\n");
        //     print_tlwe_sample(tmp, secret_keyset->tlwe_key_l1, 1./lr_params.X_beta_scale);

        //     printf(" ");
        // }
        // #endif

        LweKeySwitch(X_beta+i, ks_l1_l0, tmp);

        // #ifdef DEBUG
        // {
        //     // write_tlwe_samples("ks_out.ctxt", tmp, secret_keyset->tlwe_key_l0->params, 1);

        //     print_tlwe_sample(X_beta+i, secret_keyset->tlwe_key_l0, 1./lr_params.X_beta_scale);
        //     printf("\n");
        //     // printf("######### DEBUG END BEG #########\n");
        //     // exit(-1);
        // }
        // #endif

        del_obj(tmp);

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", i+1, lr_params.n);
    }
    VERBOSE_2_PRINT("\n");

}


/**
 * @brief Compute Xt.\sigma(X.beta)
 *
 * @param acc TRLWE sample array (n elements) where results are stored
 * @param sigmoid_xt_tps TRLWE sample array of test polynomials
 * @param trgsw_params TRGSW (and TRLWE) parameters
 * @param X_beta TLWE sample array (n elements)
 * @param tlwe_params TLWE parameters
 * @param bk_fft bootstrapping key
 */
void compute_Xt_sigma(
    TLweSample<Torus> *acc,
    const TLweSample<Torus> *sigmoid_xt_tps,
    const TGswParams<Torus>* trgsw_params,

    const LweSample<Torus>* X_beta,
    const LweParams<Torus>* tlwe_params,

    TGswSampleFFT<Torus>* bk_fft,
    const LRParams& lr_params
    )
{
    const TLweParams<Torus>* trlwe_params = trgsw_params->tlwe_params;

    /* compute vector Xt.\sigma(u).\alpha of size k (vector element are encoded in k polynomial coefficients) */
    VERBOSE_1_PRINT("==> blindrotate to compute rows of \\alpha.Xt.\\sigma(X.beta)\n");
    #pragma omp parallel for ordered schedule(static,1)
    for (int i = 0; i < lr_params.n; ++i) {
        tLweCopy(acc+i, sigmoid_xt_tps+i, trlwe_params);
        blindrotate(acc+i, trlwe_params, X_beta+i, tlwe_params, bk_fft, trgsw_params, lr_params);

        #ifdef DEBUG
        // {
        //     printf("######### DEBUG MSG BEG #########\n");
        //     printf("inp: ");
        //     print_tlwe_sample(X_beta+i, secret_keyset->tlwe_key_l0, 1./lr_params.X_beta_scale);

        //     printf("acc_%d: ", i);
        //     //double scale = 1;
        //     double scale = 1/lr_params.sigmoid_scale;
        //     print_trlwe_sample_coef(acc+i, secret_keyset->trlwe_key_l2, 0, scale);
        //     print_trlwe_sample_coef(acc+i, secret_keyset->trlwe_key_l2, 1, scale);
        //     print_trlwe_sample_coef(acc+i, secret_keyset->trlwe_key_l2, 2, scale);
        //     print_trlwe_sample_coef(acc+i, secret_keyset->trlwe_key_l2, 3, scale);
        //     printf("\n");
        //     // print_trlwe_sample(acc+i, secret_keyset->trlwe_key_l2, secret_keyset->trlwe_key_l2->params->N, scale);

        //     printf("######### DEBUG MSG END #########\n");
        // }
        #endif

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", i+1, lr_params.n);
    }
    VERBOSE_2_PRINT("\n");
}

/**
 * @brief Extracts k coefficients from \c acc and add elements from \c Xt_y_arr to them
 *  Coefficients extracted from \c acc are j for j=0..k-1
 *
 * @param beta_out TLWE sample array (k elements) to be updated
 * @param acc TRLWE sample whose k coefficients (d-spaced) are alpha.Xt.\sigma(X.beta)
 * @param Xt_y_arr TLWE sample array (k elements) containing Xt.y.alpha
 * @param trlwe_params trlwe parameters
 * @param lr_params
 */
void extract_sub_Xt_y(
    LweSample<Torus> *beta_out,
    const TLweSample<Torus> *acc,
    const LweSample<Torus>* Xt_y_arr, // k TLWE L2 encryptions of vector Xt.y.\alpha elements
    const TLweParams<Torus>* trlwe_params,
    const LRParams& lr_params
    )
{
    const LweParams<Torus>* tlwe_params = &(trlwe_params->extracted_lweparams);

    #pragma omp parallel for ordered schedule(static,1)
    for (int j = 0; j < lr_params.k; ++j) {
        LweSample<Torus>* tmp = new_obj<LweSample<Torus>>(tlwe_params);

        /* extract individual TLWE L2 samples of Xt.\sigma(u).\alpha from TRLWE L2 coefficient encoding */
        TLweFunctions<Torus>::ExtractLweSampleIndex(tmp, acc, j, tlwe_params, trlwe_params);

        /* compute delta_beta = \alpha.(Xt.\sigma(X.beta) - Xt.y)*/
        LweFunctions<Torus>::SubTo(tmp, Xt_y_arr+j, tlwe_params);

        /* update beta */
        LweFunctions<Torus>::SubTo(beta_out+j, tmp, tlwe_params);
        del_obj(tmp);

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", j+1, lr_params.n);
    }
    VERBOSE_2_PRINT("\n");
}

