#include "Solver.h"
#include<unsupported/Eigen/src/SparseExtra/MarketIO.h>
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

SparseMatrix<double> H;

void PreFE(int nelx, int nely, int* ik, int* jk)
{
    MatrixXi nodenrs(nely + 1, nelx + 1);
    int* edofVec = new int[nelx * nely];
    MatrixXi edofMat(nelx * nely, 8);
    int edofs[8] = { -1, 0, 2 * nely + 1, 2 * nely + 2, 2 * nely + 3, 2 * nely + 4, 1, 2 };

    for (size_t y = 0; y < nely + 1; y++)
    {
        for (size_t x = 0; x < nelx + 1; x++)
        {
            nodenrs(y, x) = x * (nely + 1) + y;
        }
    }

    for (size_t y = 0; y < nely; y++)
    {
        for (size_t x = 0; x < nelx; x++)
        {
            edofVec[y + x * nely] = 2 * nodenrs(y, x) + 1;
        }
    }

    for (size_t i = 0; i < nelx * nely; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            edofMat(i, j) = edofVec[i] + edofs[j];
        }
    }

    auto a = kroneckerProduct(edofMat, MatrixXi::Ones(8, 1)).eval();
    auto za = a.transpose();
    auto b = kroneckerProduct(edofMat, MatrixXi::Ones(1, 8)).eval();
    auto zb = b.transpose();

    for (size_t i = 0; i < za.cols(); i++)
    {
        for (size_t j = 0; j < za.rows(); j++)
        {
            ik[i * za.rows() + j] = za(j, i);
        }
    }

    for (size_t i = 0; i < zb.cols(); i++)
    {
        for (size_t j = 0; j < zb.rows(); j++)
        {
            jk[i * zb.rows() + j] = zb(j, i);
        }
    }

    
    delete[] edofVec;
}

void Assembly_Solve(int num_freeDofs, int num_allDofs, int num_triplets, int* free_dofs, int* ik, int* jk, double* vk, double* F, double* U)
{
    std::vector<Triplet<double>> triplets;
    triplets.reserve(num_triplets);

    for (int i = 0; i < num_triplets; i++)
    {
        triplets.push_back(Triplet<double>(ik[i], jk[i], vk[i]));
    }

    SparseMatrix<double> K(num_allDofs, num_allDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<Triplet<double>> P_triplets;
    P_triplets.reserve(num_freeDofs);

    for (int i = 0; i < num_freeDofs; i++)
    {
        P_triplets.push_back(Triplet<double>(i, free_dofs[i], 1.0));
    }

    SparseMatrix<double> P(num_freeDofs, num_allDofs);
    P.setFromTriplets(P_triplets.begin(), P_triplets.end());

    SparseMatrix<double> K_freedof = P * K * P.transpose();
    // saveMarket(K_freedof, "E:/test/K_freedof_cpp.mtx");

    VectorXd F_freedof(num_freeDofs);
    for (int i = 0; i < num_freeDofs; i++)
    {
        F_freedof(i) = F[free_dofs[i]];
    }

    VectorXd result;

    //SimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
    //llt.analyzePattern(K_freedof);
    //llt.factorize(K_freedof);
    //result = llt.solve(F_freedof);
    CholmodSimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
    llt.analyzePattern(K_freedof);
    llt.factorize(K_freedof);
    result = llt.solve(F_freedof);

    for (int i = 0; i < num_freeDofs; i++)
        U[i] = result(i);
}

double TransposeMultiply(int rows, int cols, double* A, double* U)
{
    Map<Matrix<double,Dynamic,Dynamic,RowMajor>> A_(A, rows,cols);
    Map<VectorXd> U_(U, rows);
    auto result = U_.transpose() * A_ * U_;
    return result.value();
}

void GetRowSum(int coo_length, int rows, int* ih, int* jh, double* vh, double* sh)
{
    std::vector<Triplet<double>> triplets;
    for (size_t i = 0; i < coo_length; i++)
    {
        triplets.push_back(Triplet<double>(ih[i], jh[i], vh[i]));
    }
    H.resize(rows, rows);
    H.setFromTriplets(triplets.begin(), triplets.end());
    VectorXd result = H * VectorXd::Ones(H.cols());
    for (size_t i = 0; i < result.count(); i++)
    {
        sh[i] = result[i];
    }
}

void Flt(int dc_length, double* dc, double* sh)
{
    Map<VectorXd> dc_(dc, dc_length);
    Map<VectorXd> sh_(sh, dc_length);
    VectorXd result = (H.selfadjointView<Lower>() * dc_).array() / sh_.array();
    for (size_t i = 0; i < result.count(); i++)
    {
        dc[i] = result[i];
    }
}