//==============================================================
// Copyright © Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#pragma once

#ifdef DPCPP_DLL_EXPORTS
#define DPCPP_DLL_API __declspec(dllexport)
#else
#define DPCPP_DLL_API __declspec(dllimport)
#endif

#include <vector>
#include <iostream>

typedef std::vector<int> IntVector;

DPCPP_DLL_API void InitializeVector(IntVector& a);

DPCPP_DLL_API void MultArrays(const long long* a, const long long* b, long long* res, int num, int nc);
