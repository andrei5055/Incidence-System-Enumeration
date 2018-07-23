#include <iostream>
#include "stdafx.h"
#include "matrix.h"

template class C_tDesign<MATRIX_ELEMENT_TYPE>;

template<class T>
C_tDesign<T>::C_tDesign(int t, int v, int k, int lambda) : C_BIBD(v, k, t), m_t(t)
{
	// Define all lambdas: 
	int i = t;
	int *pLambda = new int[t - 2];
	while (--i > 1) {
		pLambda[i - 2] = lambda;
		lambda = lambda * (v - i) / (k - i);
	}

	// Initiate BIBD's parameter
	Init_BIBD_param(v, k, lambda);

	// Add remaining lambda's to Lambda set
	while (++i < t)
		AddValueToNumSet(pLambda[i - 2], t_lSet);

	// Initiate matrix 
	Init(v, lambda * v * (v - 1) / (k * (k - 1)));
	delete[] pLambda;
}
