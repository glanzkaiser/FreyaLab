#pragma once
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>

// High-performance Real Compressed Sparse Row (CSR) Storage
struct CsrMatrix {
    int numRows;
    std::vector<double> values;
    std::vector<int> colIndices;
    std::vector<int> rowPointers; // Size: numRows + 1

    void initializeFromPattern(const CsrMatrix& referenceMatrix); 
    void initialize(int rows, const std::vector<std::pair<int, int>>& sparsityPattern);
    void zero();
    void addValue(int row, int col, double val);
};

// High-performance Complex Compressed Sparse Row (CSR) Storage
struct CsrMatrixComplex {
    int numRows;
    std::vector<std::complex<double>> values;
    std::vector<int> colIndices;
    std::vector<int> rowPointers; // Size: numRows + 1

    void initializeFromPattern(const CsrMatrix& referenceMatrix);
    void zero();
    void addValue(int row, int col, const std::complex<double>& val);
};

// Advanced Numerical Linear Algebra & SPICE Limiting Solvers
class CircuitMath {
public:
    // --- Real-Valued Left-Looking Gilbert-Peierls LU Factorization Architecture ---
    static bool luFactorize(const CsrMatrix& A, CsrMatrix& L, CsrMatrix& U, std::vector<int>& P);
    static std::vector<double> luSolve(const CsrMatrix& L, const CsrMatrix& U, const std::vector<int>& P, const std::vector<double>& B);
    static void gpDfs(int j, const std::vector<std::vector<int>>& u_rows_cols, 
                      std::vector<bool>& visited, std::vector<int>& xi, int& top, int i);

    // --- Complex-Valued Left-Looking Gilbert-Peierls LU Factorization Architecture ---
    static bool luFactorizeComplex(const CsrMatrixComplex& A, CsrMatrixComplex& L, CsrMatrixComplex& U, std::vector<int>& P);
    static std::vector<std::complex<double>> luSolveComplex(const CsrMatrixComplex& L, const CsrMatrixComplex& U, const std::vector<int>& P, const std::vector<std::complex<double>>& B);
    static void gpDfsComplex(int j, const std::vector<std::vector<int>>& u_rows_cols, 
                         std::vector<bool>& visited, std::vector<int>& xi, int& top, int i);


    // --- SPICE Non-linear Voltage Limiting (Modified PN Junction Limiting) ---
    static double pnjlim(double vnew, double vold, double vt, double vcrit);
};

