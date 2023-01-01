#include "Solver.h"
#include<unsupported/Eigen/src/SparseExtra/MarketIO.h>
#include <unsupported/Eigen/KroneckerProduct>
#include<ctime>

using namespace Eigen;
using namespace std;

SparseMatrix<double> H;

void PrintMatrix(int row, int col, double* val)
{
    MatrixXd K(row, col);
    for (size_t i = 0; i < row; i++)
    {
        for (size_t j = 0; j < col; j++)
        {
            K(i, j) = val[i * col + j];
        }
    }
    cout << K << std::endl;
}


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

void Assembly_Solve(int num_freeDofs, int num_allDofs, int num_triplets, int* free_dofs, int* ik, int* jk, double* sk, double* F, double* U)
{
    std::vector<Triplet<double>> triplets;
    triplets.reserve(num_triplets);

    std::cout << sk[1500] << std::endl;
    for (int i = 0; i < num_triplets; i++)
    {
        triplets.push_back(Triplet<double>(ik[i], jk[i], sk[i]));
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

    auto KD = MatrixXd(K_freedof);
    std::cout << KD(94, 94) << std::endl;
    //Eigen::saveMarket(K_freedof, "mat.mtx");

    VectorXd F_freedof(num_freeDofs);
    for (int i = 0; i < num_freeDofs; i++)
    {
        F_freedof(i) = F[free_dofs[i]];
    }

    VectorXd result;

    clock_t _start;
    clock_t _end;
    _start = clock();
    PardisoLLT<Eigen::SparseMatrix<double>,1> llt(K_freedof);
    //llt.pardisoParameterArray()[59] = 0;
    //mkl_set_num_threads(1);
    //CholmodSimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
    //CholmodSupernodalLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
    llt.analyzePattern(K_freedof);
    llt.factorize(K_freedof);
    result = llt.solve(F_freedof);
    _end = clock();
    double endtime = (double)(_end - _start) / CLOCKS_PER_SEC;
    cout << "Solving cost:" << endtime * 1000 << "ms" << endl;	//msÎªµ¥Î»
    Eigen::VectorXd::Map(U, result.rows()) = result;
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

    Eigen::VectorXd::Map(sh, result.rows()) = result;
}

void Flt(int dc_length, double* dc, double* sh)
{
    Map<VectorXd> dc_(dc, dc_length);
    Map<VectorXd> sh_(sh, dc_length);
    VectorXd result = (H.selfadjointView<Lower>() * dc_).array() / sh_.array();    
    Eigen::VectorXd::Map(dc, result.rows()) = result;
}