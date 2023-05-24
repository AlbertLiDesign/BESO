#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <chrono>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Eigen>
#include <Eigen/PardisoSupport>

#include <unsupported/Eigen/SparseExtra>
#include<unsupported/Eigen/src/SparseExtra/MarketIO.h>
#include <unsupported/Eigen/KroneckerProduct>

#include<ctime>
#include<time.h>


extern "C" __declspec(dllexport) void Assembly_Solve(bool parallel, int num_freeDofs, int num_allDofs, int num_triplets,
	int* free_dofs, int* ik, int* jk, double* sk, double* F, double* U);
extern "C" __declspec(dllexport) void PreFE(int nelx, int nely, int* ik, int* jk);
extern "C" __declspec(dllexport) double TransposeMultiply(int rows, int cols, double* A, double* U);
extern "C" __declspec(dllexport) void GetRowSum(int coo_length, int rows, int* ih, int* jh, double* vh, double* sh);
extern "C" __declspec(dllexport) void Flt(int dc_length, double* dc, double* sh);
extern "C" __declspec(dllexport) void Flt3D(int nEl, double* dc, double* sh);
extern "C" __declspec(dllexport) void PrintMatrix(int row, int col, double* val);
extern "C" __declspec(dllexport) void Cal_ik_jk(int nEl, int* edofMat, int* ik, int* jk);