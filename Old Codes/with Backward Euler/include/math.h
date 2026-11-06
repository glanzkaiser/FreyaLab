#pragma once

#include <vector>
#include <iostream>
#include <utility> 

struct CsrMatrix {
    int numRows = 0;
    std::vector<double> values;     
    std::vector<int> colIndices;    
    std::vector<int> rowPointers;   

    void clearNumericalValues() {
        std::fill(values.begin(), values.end(), 0.0);
    }

    void clear() {
        values.clear();
        colIndices.clear();
        rowPointers.clear();
        numRows = 0;
    }
    
    void reset() {
        clear();
    }
};

class MathSolver {
public:
    static void accumulateCsr(CsrMatrix& A, int row, int col, double val);
    static void finalizeStructure(CsrMatrix& A, const std::vector<std::pair<int, int>>& nonZeroPatterns);
    static bool sparseLuFactorize(const CsrMatrix& A, CsrMatrix& L, CsrMatrix& U, std::vector<int>& pivot);
    static void sparseLuSolve(const CsrMatrix& L, const CsrMatrix& U, const std::vector<int>& pivot,
                              const std::vector<double>& B, std::vector<double>& x);
    static double pnjlim(double vnew, double vold, double vt, double vcrit);
};



