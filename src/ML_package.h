#ifdef MATLAB_MEX_FILE

#ifndef ML_PACKAGE_H
#define ML_PACKAGE_H

#ifdef ES_DEV_ENV
	#include "../../../SMat/smat.def.h"
	#include "../../../SMat/smat.h"
	#include "../../../MLDev/ML_meal.h"
#else
	#include "smat.def.h"
	#include "smat.h"
	#include "ML_meal.h"
#endif

__declspec(dllexport) void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[]);


void ML_tkmeans (int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[])

#endif	//	#ifndef ML_PACKAGE_H
#endif	//	#ifdef MATLAB_MEX_FILE
