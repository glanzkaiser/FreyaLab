// freyalab/KLU.h
#pragma once

#include "freyalab/math.h"
#include <complex>
#include <vector>

struct KluBlockFactor {
    int startIdx = 0;
    int size = 0;

    // Local fill-reducing AMD permutation inside this BTF block
    std::vector<int> amdOrder;      // local new index -> local BTF index
    std::vector<int> amdInvOrder;   // local BTF index -> local new index

    // Partial pivoting & pre-allocated CSR L and U factors for real-valued solves
    std::vector<int> localPivotRow;
    CsrMatrix L;                    // Unit lower triangular (strictly below diagonal)
    CsrMatrix U;                    // Upper triangular (strictly above diagonal)
    std::vector<double> diagU;      // Cached real U diagonal for O(1) division

    // Pre-allocated CSR L and U factors for complex-valued AC phasor solves
    CsrMatrixComplex L_c;
    CsrMatrixComplex U_c;
    std::vector<std::complex<double>> diagU_c;
};

struct KluSymbolicContext {
    int n = 0;
    bool isAnalyzed = false;

    // Global BTF row/column permutations
    std::vector<int> rowPerm;       // BTF row index -> original matrix row
    std::vector<int> invRowPerm;    // original matrix row -> BTF row index
    std::vector<int> colPerm;       // BTF col index -> original matrix col
    std::vector<int> invColPerm;    // original matrix col -> BTF col index

    // Diagonal BTF blocks (solved in reverse order during back-substitution)
    std::vector<KluBlockFactor> blocks;

    // Off-diagonal upper-block CSR matrices (in BTF coordinates) for block back-substitution
    CsrMatrix offDiagUpper;
    CsrMatrixComplex offDiagUpper_c;

    // Pre-allocated zero-allocation real numeric workspaces
    int maxBlockSize = 0;
    std::vector<double> workX;
    std::vector<int> visitedEpoch;
    std::vector<int> touchedEpoch;
    std::vector<int> touchedCols;
    std::vector<int> xi;
    int currentEpoch = 1;

    // Reusable real solve scratch buffers
    std::vector<double> btfRhs;
    std::vector<double> btfSol;

    // Pre-allocated complex numeric workspaces & scratch buffers
    std::vector<std::complex<double>> workX_c;
    std::vector<std::complex<double>> btfRhs_c;
    std::vector<std::complex<double>> btfSol_c;
};

class KluSolver {
public:
    // Threshold where KLU (BTF + AMD + Symbolic Cache) outperforms direct CSR LU
    static constexpr int KLU_ACTIVATION_THRESHOLD = 300;

    // One-time symbolic analysis: Maximum Transversal + Tarjan BTF + AMD Ordering
    static KluSymbolicContext analyzeSymbolic(const CsrMatrix& pattern);

    // Fast zero-allocation real numeric factorization reusing pre-analyzed BTF + AMD structure
    static bool factorizeNumeric(const CsrMatrix& A, KluSymbolicContext& sym);

    // Fast real block forward/backward triangular solve: A * x = b
    static std::vector<double> solveNumeric(KluSymbolicContext& sym, const std::vector<double>& b);

    // Fast complex numeric factorization reusing the same pre-analyzed BTF + AMD structure
    static bool factorizeNumericComplex(const CsrMatrixComplex& A, KluSymbolicContext& sym);

    // Fast complex block forward/backward triangular solve: A_c * x_c = b_c
    static std::vector<std::complex<double>> solveNumericComplex(KluSymbolicContext& sym,
                                                                 const std::vector<std::complex<double>>& b);
};

