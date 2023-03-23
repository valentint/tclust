//  VT::27.06.2017 - this file was added to fix the warning "Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'"
//
//  About registration of native symbols see for example: https://www.r-bloggers.com/1-easy-package-registration/
//      also here http://r.789695.n4.nabble.com/Registration-of-native-routines-td4728874.html
//      - about Windows - take the 64 bit version of mingw!
//

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include <R.h>              
#include <Rinternals.h>

/* .C calls */
/*
EXPORT void RestrictEigenValues			(int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize) ;
EXPORT void RestrictEigenValues_deter	(int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdEv, double *pdClustSize) ;
EXPORT void C_tclust (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, double *pdS, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdZ, double *pdObjER, int *pnConvER) ;
EXPORT void C_tkmeans (int *pnParIn, int *pnParOut, double *pdParIn, double *pdParOut, double *pdX, double *pdM, int *pnAssign, double *pdClustSize, double *pdWeights, double *pdObjER, int *pnConvER) ;

EXPORT void restr_dir_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval) ;
EXPORT void restr_prop_C (int *pdwParamIn, int *pdwParamOut, double *pdParamIn, double *pdS, double *pdCSize, double *pdEVec, double *pdEval) ;
EXPORT void OptVectors_C (int *pdwParamIn, double *pdU, double *pdD, double *pdS, double *pdCSize) ;
EXPORT void sme_eigen_sym_2x2_norm (double * const pdEval, double *const pdEVec, const double *const pd, const double *pdZeroTol) ;
*/

extern void RestrictEigenValues(void *, void *, void *, void *, void *, void *);
extern void RestrictEigenValues_deter(void *, void *, void *, void *, void *, void *);
extern void C_tclust(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_tkmeans(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static R_NativePrimitiveArgType RestrictEigenValues_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType RestrictEigenValues_deter_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_tclust_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType C_tkmeans_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

static const R_CMethodDef CEntries[] = {
    {"RestrictEigenValues",       (DL_FUNC) &RestrictEigenValues,        6, RestrictEigenValues_t},
    {"RestrictEigenValues_deter", (DL_FUNC) &RestrictEigenValues_deter,  6, RestrictEigenValues_deter_t},
    {"C_tclust",                    (DL_FUNC) &C_tclust,                13, C_tclust_t},
    {"C_tkmeans",                   (DL_FUNC) &C_tkmeans,               11, C_tkmeans_t},
    {NULL, NULL, 0}
};

void R_init_tclust(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE); 
}
