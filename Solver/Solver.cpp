#include "Solver.h"

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

void Cal_ik_jk(int nEl, int* edofMat, int* ik, int* jk)
{
    MatrixXi A = Map<MatrixXi>(edofMat, nEl, 24);

    auto a = kroneckerProduct(A, MatrixXi::Ones(24, 1)).eval();
    auto za = a.transpose();
    auto b = kroneckerProduct(A, MatrixXi::Ones(1, 24)).eval();
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
}
void PreFE(int nelx, int nely, int* ik, int* jk)
{
    MatrixXi nodenrs(nely + 1, nelx + 1);
    int* edofVec = new int[nelx * nely];
    MatrixXi edofMat(nelx * nely, 8);
    int edofs[8] = { -1, 0, 2 * nely + 1, 2 * nely + 2, 2 * nely + 3, 2 * nely + 4, 1, 2 };

    for (int y = 0; y < nely + 1; y++)
    {
        for (int x = 0; x < nelx + 1; x++)
        {
            nodenrs(y, x) = x * (nely + 1) + y;
        }
    }
    for (int y = 0; y < nely; y++)
    {
        for (int x = 0; x < nelx; x++)
        {
            edofVec[y + x * nely] = 2 * nodenrs(y, x) + 1;
        }
    }
    for (int i = 0; i < nelx * nely; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            edofMat(i, j) = edofVec[i] + edofs[j];
        }
    }
    auto a = kroneckerProduct(edofMat, MatrixXi::Ones(8, 1)).eval();
    auto za = a.transpose();
    auto b = kroneckerProduct(edofMat, MatrixXi::Ones(1, 8)).eval();
    auto zb = b.transpose();
    for (int i = 0; i < za.cols(); i++)
    {
        for (int j = 0; j < za.rows(); j++)
        {
            ik[i * za.rows() + j] = za(j, i);
        }
    }
    for (int i = 0; i < zb.cols(); i++)
    {
        for (int j = 0; j < zb.rows(); j++)
        {
            jk[i * zb.rows() + j] = zb(j, i);
        }
    }
    delete[] edofVec;
}

void Assembly_Solve(int solver, bool parallel, int num_freeDofs, int num_allDofs, int num_triplets, int* free_dofs, int* ik, int* jk, double* sk, double* F, double* U)
{
    //auto start2 = std::chrono::high_resolution_clock::now();
    std::vector<Triplet<double>> triplets(num_triplets);

    for (int i = 0; i < num_triplets; i++)
    {
        triplets[i] = Triplet<double>(ik[i], jk[i], sk[i]);
    }

    SparseMatrix<double> K(num_allDofs, num_allDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<Triplet<double>> P_triplets(num_freeDofs);

    for (int i = 0; i < num_freeDofs; i++)
    {
        P_triplets[i] = Triplet<double>(i, free_dofs[i], 1.0);
    }

    SparseMatrix<double> P(num_freeDofs, num_allDofs);
    P.setFromTriplets(P_triplets.begin(), P_triplets.end());
    auto end2 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> elapsed2 = end2 - start2;
    //std::cout << "Assembly: " << elapsed2.count() << std::endl;

    SparseMatrix<double> K_freedof = P * K * P.transpose();
    //Eigen::saveMarket(K_freedof, "K_freedof.mtx");

    VectorXd F_freedof(num_freeDofs);
    for (int i = 0; i < num_freeDofs; i++)
        F_freedof(i) = F[free_dofs[i]];


    VectorXd result;
    PardisoLLT<Eigen::SparseMatrix<double>, 1> pardiso;
    SimplicialLLT<Eigen::SparseMatrix<double>> simplicial;
    CholmodSupernodalLLT<Eigen::SparseMatrix<double>> supernodal;
    switch (solver)
    {
    case 0:
        simplicial.analyzePattern(K_freedof); // 预先分析A的非零模式
        simplicial.factorize(K_freedof); // 使用已分析的模式进行分解
        result = simplicial.solve(F_freedof); // 使用分解后的矩阵求解)
        break;
    case 1:
        if (!parallel)
        {
            pardiso.pardisoParameterArray()[59] = 0;
            mkl_set_num_threads(1);
        }

        pardiso.analyzePattern(K_freedof);
        pardiso.factorize(K_freedof);
        result = pardiso.solve(F_freedof);
        break;
    case 2:
        supernodal.analyzePattern(K_freedof); // 预先分析A的非零模式
        supernodal.factorize(K_freedof); // 使用已分析的模式进行分解
        result = supernodal.solve(F_freedof); // 使用分解后的矩阵求解
        break;
    default:
        simplicial.analyzePattern(K_freedof); // 预先分析A的非零模式
        simplicial.factorize(K_freedof); // 使用已分析的模式进行分解
        result = simplicial.solve(F_freedof); // 使用分解后的矩阵求解)
        break;
    }   

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

void Flt3D(int nEl, double* dc, double* sh)
{
    VectorXd A(nEl);
    for (size_t i = 0; i < nEl; i++)
        A(i) = dc[i] / sh[i];
    VectorXd result = H.selfadjointView<Lower>() * A;
    Eigen::VectorXd::Map(dc, result.rows()) = result;
}