#include "lr_params.h"
#include "keyset.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"


/**
 * @brief Blind rotate \c acc by TLWE sample \c input coefficients
 */
void blindrotate(
    TLweSample<Torus>* acc,
    const TLweParams<Torus>* acc_params,
    const LweSample<Torus>* inp,
    const LweParams<Torus>* inp_params,
    const TGswSampleFFT<Torus>* bk_fft,
    const TGswParams<Torus>* bk_params
    )
{
    const int N=acc_params->N;
    const int _2N = 2*N;
    const int n = inp_params->n;

    // Modulus switching
    int* bara = new int[N];
    int barb = TorusUtils<Torus>::modSwitchFromTorus(inp->b,_2N);
    for (int i=0; i<n; i++) {
        bara[i]=TorusUtils<Torus>::modSwitchFromTorus(inp->a[i],_2N);
    }

    TLweSample<Torus>* tmp = new_obj<TLweSample<Torus>>(acc_params);
    tLweCopy(tmp, acc, acc_params);

    // Accumulator rotate by -barb
    TLweFunctions<Torus>::MulByXai(acc, _2N-barb, tmp, acc_params);

    // Accumulator blind rotate by \sum bara_i * sk_i
    tfhe_blindRotate_FFT(acc, bk_fft, bara, n, bk_params);

    delete[] bara;
}



/**
 * @brief Compute Xt.y from TRGSW samples of X columns and TRLWE sample of y (coefficient packed)
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
        TLweFunctions<Torus>::ExtractLweSampleIndex(Xt_y+i, acc, lr_params.n, tlwe_params, trlwe_params);
    }
}


/**
 * @brief Compute Xt.(1/2 - y)
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
    TorusPolynomial<Torus>* yp_poly = new_obj<TorusPolynomial<Torus>>(trlwe_params->N);
    Torus mu = TorusUtils<Torus>::from_double(1./2);
    for (int i = 0; i < lr_params.n; ++i)
        yp_poly->coefsT[i] = mu;

    TLweSample<Torus>* yp = new_obj<TLweSample<Torus>>(trlwe_params);
    TLweFunctions<Torus>::NoiselessTrivial(yp, yp_poly, trlwe_params);
    del_obj(yp_poly);

    TLweFunctions<Torus>::SubTo(yp, y, trlwe_params);

    compute_Xt_y(beta, X_cols, yp, trgsw_params, lr_params);
    del_obj(yp);
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
    const TGswParams<Torus>* trgsw_params,
    const LRParams& lr_params
    )
{
    /* switch beta from TLWE to TRLWE and multiply each beta by X columns */
    VERBOSE_1_PRINT("==> compute individual columns of X.beta\n");
    #pragma omp parallel for ordered schedule(static,1)
    for (int j = 0; j < lr_params.k; ++j) {
        /* switch from L2 TLWE sample to L1 TRLWE sample  */
        //KeySwitch(X_beta_cp+j, keyset->ks_l2_l1, beta_in+j);

        TGswFunctions<Torus>::ExternMulToTLwe(X_beta_cp+j, X_cols+j, trgsw_params);

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
        LweFunctions<Torus>::KeySwitch(X_beta+i, ks_l1_l0, tmp);
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
        blindrotate(acc+i, trlwe_params, X_beta+i, tlwe_params, bk_fft, trgsw_params);

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", i+1, lr_params.n);
    }
    VERBOSE_2_PRINT("\n");
}
/**
 * @brief Accumulate array elements and store result to first element
 *
 * @param acc TRLWE sample array (n elements)
 * @param trlwe_params TRLWE sample parameters
 * @param lr_params
 */
void accumulate_trlwe_acc(
    TLweSample<Torus> *acc,
    const TLweParams<Torus>* trlwe_params,
    const LRParams& lr_params
    )
{
    VERBOSE_1_PRINT("==> accumulate to compute \\alpha.Xt.\\sigma(X.beta) ... ");
    for (int i = 1; i < lr_params.n; ++i) {
        TLweFunctions<Torus>::AddTo(acc+0, acc+i, trlwe_params);
    }
    VERBOSE_1_PRINT("done\n");
}

/**
 * @brief Extracts k coefficients from \c acc and add elements from \c Xt_y_arr to them
 *  Coefficients extracted from \c acc are j*lr_params.d for j=0..k-1
 *
 * @param beta_out TLWE sample array (k elements) to be updated
 * @param acc TRLWE sample whose k coefficients (d-spaced) are alpha.Xt.\sigma(X.beta)
 * @param Xt_y_arr TLWE sample array (k elements) containing Xt.y.alpha
 * @param trlwe_params trlwe parameters
 * @param lr_params
 */
void extract_add_Xt_y(
    LweSample<Torus> *beta_out,
    const TLweSample<Torus> *acc,
    const LweSample<Torus>* Xt_y_arr, // k TLWE L2 encryptions of vector Xt.y.\alpha elements
    const TLweParams<Torus>* trlwe_params,
    const LRParams& lr_params
    )
{
    const LweParams<Torus>* tlwe_params = &(trlwe_params->extracted_lweparams);

    VERBOSE_1_PRINT("==> compute \\alpha.Xt.[\\sigma(X.beta)-y]\n");
    #pragma omp parallel for ordered schedule(static,1)
    for (int j = 0; j < lr_params.k; ++j) {
        LweSample<Torus>* tmp = new_obj<LweSample<Torus>>(tlwe_params);

        /* extract individual TLWE L2 samples of Xt.\sigma(u).\alpha from TRLWE L2 coefficient encoding */
        TLweFunctions<Torus>::ExtractLweSampleIndex(tmp, acc, j*lr_params.d, tlwe_params, trlwe_params);
        /* compute delta_beta = \alpha.(Xt.\sigma(X.beta) - Xt.y)*/
        LweFunctions<Torus>::SubTo(tmp, Xt_y_arr+j, tlwe_params);

        /* update beta */
        LweFunctions<Torus>::AddTo(beta_out+j, tmp, tlwe_params);
        del_obj(tmp);

        #pragma omp ordered
        VERBOSE_2_PRINT("     %4d/%4d\r", j+1, lr_params.n);
    }
    VERBOSE_2_PRINT("\n");
}
