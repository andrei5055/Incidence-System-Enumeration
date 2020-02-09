#include <iostream>
#include <algorithm>
#include "stdafx.h"
#include "matrix.h"

template class C_tDesign<MATRIX_ELEMENT_TYPE>;
template class CCombinedBIBD<MATRIX_ELEMENT_TYPE>;

template<class T>
C_tDesign<T>::C_tDesign(int t, int v, int k, int lambda) : C_BIBD<T>(v, k, t), m_t(t)
{
	// Define all lambdas: 
	int i = t;
	int *pLambda = new int[t - 2];
	while (--i > 1) {
		pLambda[i - 2] = lambda;
		lambda = lambda * (v - i) / (k - i);
	}

	// Initiate BIBD's parameter
	this->InitParam(v, k, lambda);

	// Add remaining lambda's to Lambda set
	while (++i < t)
		this->AddValueToNumSet(pLambda[i - 2], t_lSet);

	// Initiate matrix 
	this->Init(v, lambda * v * (v - 1) / (k * (k - 1)));
	delete[] pLambda;
}


template<class T>
CCombinedBIBD<T>::CCombinedBIBD(int v, int k, const std::vector<int>& lambdaInp) : C_BIBD<T>(0, k) {
	std::vector<int> lambdaSet(lambdaInp);
	std::sort(lambdaSet.begin(), lambdaSet.end(), std::greater<int>());
	m_ppParamSet = createParamStorage(t_rSet); // Create 2 sets of vector (for Lambda and R of the component of combineed BIBD)

	const auto nSubDesigns = lambdaSet.size();
	const auto v1 = v - 1;
	const auto k1 = k - 1;
	int lambda = 0;
	for (size_t i = 0; i < nSubDesigns; ++i) {
		const auto lambdaCurr = lambdaSet[i];
		const auto r = lambdaCurr * v1 / k1;
		assert(r * k1 == lambdaCurr * v1);
		m_ppParamSet[t_lSet]->AddElement(lambdaCurr);
		m_ppParamSet[t_rSet]->AddElement(r);
		lambda += lambdaCurr;
	}

	Init(v + 1, lambda * v * v1 / (k * k1));
	InitParam(v, k, lambda);
}