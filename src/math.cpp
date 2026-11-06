#include "freyalab/math.h"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>

// ============================================================================
// REAL CSR MATRIX STORAGE ENGINE
// ============================================================================
//good and robust
void CsrMatrix::initialize(int rows, const std::vector<std::pair<int, int>>& sparsityPattern) {
    numRows = rows;
    std::vector<std::set<int>> sortedPattern(rows);
    for (const auto& entry : sparsityPattern) {
        if (entry.first >= 0 && entry.first < rows && entry.second >= 0 && entry.second < rows) {
            sortedPattern[entry.first].insert(entry.second);
        }
    }

    rowPointers.clear();
    colIndices.clear();
    values.clear();

    rowPointers.push_back(0);
    for (int i = 0; i < rows; ++i) {
        for (int col : sortedPattern[i]) {
            colIndices.push_back(col);
            values.push_back(0.0);
        }
        rowPointers.push_back(static_cast<int>(colIndices.size()));
    }
}
void CsrMatrix::initializeFromPattern(const CsrMatrix& referenceMatrix) {
    numRows = referenceMatrix.numRows;
    rowPointers = referenceMatrix.rowPointers;
    colIndices = referenceMatrix.colIndices;
    values.assign(colIndices.size(), 0.0);
}


void CsrMatrix::zero() {
    std::fill(values.begin(), values.end(), 0.0);
}

void CsrMatrix::addValue(int row, int col, double val) {
    if (row < 0 || row >= numRows) return;
    int start = rowPointers[row];
    int end = rowPointers[row + 1];
    for (int i = start; i < end; ++i) {
        if (colIndices[i] == col) {
            values[i] += val;
            return;
        }
    }
}

// ============================================================================
// COMPLEX CSR MATRIX STORAGE ENGINE
// ============================================================================
//good and robust
void CsrMatrixComplex::initializeFromPattern(const CsrMatrix& referenceMatrix) {
    numRows = referenceMatrix.numRows;
    rowPointers = referenceMatrix.rowPointers;
    colIndices = referenceMatrix.colIndices;
    values.assign(colIndices.size(), std::complex<double>(0.0, 0.0));
}

void CsrMatrixComplex::zero() {
    std::fill(values.begin(), values.end(), std::complex<double>(0.0, 0.0));
}

void CsrMatrixComplex::addValue(int row, int col, const std::complex<double>& val) {
    if (row < 0 || row >= numRows) return;
    int start = rowPointers[row];
    int end = rowPointers[row + 1];
    for (int i = start; i < end; ++i) {
        if (colIndices[i] == col) {
            values[i] += val;
            return;
        }
    }
}

// ============================================================================
// REAL-VALUED LINEAR ALGEBRA TRACK
// ============================================================================
// Fix: Updated signature to track the upper-triangular row dependency matrix 'U' up to row index 'i'
void CircuitMath::gpDfs(int j, const std::vector<std::vector<int>>& u_rows_cols, std::vector<bool>& visited, std::vector<int>& xi, int& top, int i) {
    visited[j] = true;
    for (int col : u_rows_cols[j]) {
        // Fix: Add explicit column constraint filtering to isolate previous topological rows
        if (col < i && !visited[col]) {
            gpDfs(col, u_rows_cols, visited, xi, top, i);
        }
    }
    xi[--top] = j;
}
// good and robust
bool CircuitMath::luFactorize(const CsrMatrix& A, CsrMatrix& L, CsrMatrix& U, std::vector<int>& P) {
    int n = A.numRows;
    P.resize(n);
    for (int i = 0; i < n; ++i) P[i] = i;

    std::vector<double> x(n, 0.0);
    std::vector<bool> visited(n, false);
    std::vector<int> xi(n);

    std::vector<std::vector<int>> l_rows_cols(n);
    std::vector<std::vector<double>> l_rows_vals(n);
    std::vector<std::vector<int>> u_rows_cols(n);
    std::vector<std::vector<double>> u_rows_vals(n);
    
    std::vector<double> current_pivots(n, 0.0);

    for (int i = 0; i < n; ++i) {
        int top = n;
        std::fill(visited.begin(), visited.end(), false);

        int start = A.rowPointers[i];
        int end = A.rowPointers[i + 1];
        for (int idx = start; idx < end; ++idx) {
            int col = A.colIndices[idx];
            x[col] = A.values[idx];
            
            if (col < i && !visited[col]) {
                gpDfs(col, u_rows_cols, visited, xi, top, i);
            }
        }

        // Numeric elimination phase using the topological ordering found
        for (int idx = top; idx < n; ++idx) {
            int j = xi[idx];
            double pivot = current_pivots[j];
            
            if (std::abs(pivot) < 1e-12) {
                pivot = (pivot < 0.0) ? -1e-12 : 1e-12;
            }

            double multiplier = x[j] / pivot;
            
            // FIX: Immediately commit the multiplier into L's sparse tracking vectors.
            // This isolates it from being overwritten by subsequent elimination updates to x[j].
            l_rows_cols[i].push_back(j);
            l_rows_vals[i].push_back(multiplier);

            for (size_t k_idx = 0; k_idx < u_rows_cols[j].size(); ++k_idx) {
                int col = u_rows_cols[j][k_idx];
                if (col > j) { 
                    x[col] -= multiplier * u_rows_vals[j][k_idx];
                }
            }
        }

        if (std::abs(x[i]) < 1e-12) {
            x[i] = (x[i] < 0.0) ? -1e-12 : 1e-12;
        }
        
        current_pivots[i] = x[i];

        // Extract and permanently store row i's structural upper-triangular components
        for (int col = i; col < n; ++col) {
            if (std::abs(x[col]) > 1e-15 || col == i) {
                u_rows_cols[i].push_back(col);
                u_rows_vals[i].push_back(x[col]);
            }
            x[col] = 0.0;
        }
        
        // Clear remaining lower-triangular working garbage elements in x
        for (int idx = top; idx < n; ++idx) {
            x[xi[idx]] = 0.0;
        }
    }

    // Assemble final CSR data matrices for L and U
    L.numRows = n; L.rowPointers.assign(n + 1, 0); L.colIndices.clear(); L.values.clear();
    U.numRows = n; U.rowPointers.assign(n + 1, 0); U.colIndices.clear(); U.values.clear();

    for (int i = 0; i < n; ++i) {
        // Sort L's row elements by column index to guarantee strict CSR compliance
        std::vector<std::pair<int, double>> l_row_elements;
        for (size_t j = 0; j < l_rows_cols[i].size(); ++j) {
            l_row_elements.push_back({l_rows_cols[i][j], l_rows_vals[i][j]});
        }
        std::sort(l_row_elements.begin(), l_row_elements.end());

        for (const auto& elem : l_row_elements) {
            L.colIndices.push_back(elem.first);
            L.values.push_back(elem.second);
        }
        L.colIndices.push_back(i);
        L.values.push_back(1.0); // Unit lower triangular matrix diagonal entry
        L.rowPointers[i + 1] = static_cast<int>(L.colIndices.size());

        for (size_t j = 0; j < u_rows_cols[i].size(); ++j) {
            U.colIndices.push_back(u_rows_cols[i][j]);
            U.values.push_back(u_rows_vals[i][j]);
        }
        U.rowPointers[i + 1] = static_cast<int>(U.colIndices.size());
    }

    return true;
}
//good and robust
std::vector<double> CircuitMath::luSolve(const CsrMatrix& L, const CsrMatrix& U, const std::vector<int>& P, const std::vector<double>& B) {
    int n = L.numRows;
    std::vector<double> x(n, 0.0), y(n, 0.0);

    // 1. Forward Substitution: L * y = P * B
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        int start = L.rowPointers[i];
        int end = L.rowPointers[i + 1]; // Fix: Process the entire row slice safely
        
        for (int j = start; j < end; ++j) {
            int col = L.colIndices[j];
            if (col == i) {
                continue; // Fix: Filter out L's diagonal elements cleanly (since L is unit lower triangular)
            }
            sum += L.values[j] * y[col];
        }
        y[i] = B[P[i]] - sum;
    }

    // 2. Backward Substitution: U * x = y
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        int start = U.rowPointers[i];
        int end = U.rowPointers[i + 1];
        double diagVal = 1.0;
        
        for (int j = start; j < end; ++j) {
            if (U.colIndices[j] == i) {
                diagVal = U.values[j];
            } else {
                sum += U.values[j] * x[U.colIndices[j]];
            }
        }
        
        if (std::abs(diagVal) < 1e-12) {
            diagVal = (diagVal < 0.0) ? -1e-12 : 1e-12; // Protect against division-by-zero
        }
        x[i] = (y[i] - sum) / diagVal;
    }

    return x;
}

// ============================================================================
// COMPLEX-VALUED LINEAR ALGEBRA TRACK
// ============================================================================
// Fix: Updated signature to look up the upper-triangular structure 'u_rows_cols' dynamically
void CircuitMath::gpDfsComplex(int j, const std::vector<std::vector<int>>& u_rows_cols, std::vector<bool>& visited, std::vector<int>& xi, int& top, int i) {
    visited[j] = true;
    for (int col : u_rows_cols[j]) {
        if (col < i && !visited[col]) {
            gpDfsComplex(col, u_rows_cols, visited, xi, top, i);
        }
    }
    xi[--top] = j;
}
//good and robust
bool CircuitMath::luFactorizeComplex(const CsrMatrixComplex& A, CsrMatrixComplex& L, CsrMatrixComplex& U, std::vector<int>& P) {
    int n = A.numRows;
    P.resize(n);
    for (int i = 0; i < n; ++i) P[i] = i;

    std::vector<std::complex<double>> x(n, std::complex<double>(0.0, 0.0));
    std::vector<bool> visited(n, false);
    std::vector<int> xi(n);

    std::vector<std::vector<int>> l_rows_cols(n);
    std::vector<std::vector<std::complex<double>>> l_rows_vals(n);
    std::vector<std::vector<int>> u_rows_cols(n);
    std::vector<std::vector<std::complex<double>>> u_rows_vals(n);
    
    // FIX 1: Dense tracker to isolate complex pivot calculation states securely
    std::vector<std::complex<double>> current_pivots(n, std::complex<double>(0.0, 0.0));

    for (int i = 0; i < n; ++i) {
        int top = n;
        std::fill(visited.begin(), visited.end(), false);

        int start = A.rowPointers[i];
        int end = A.rowPointers[i + 1];
        for (int idx = start; idx < end; ++idx) {
            int col = A.colIndices[idx];
            x[col] = A.values[idx];
            
            if (col < i && !visited[col]) {
                gpDfsComplex(col, u_rows_cols, visited, xi, top, i);
            }
        }

        // Numeric elimination phase using the topological ordering found
        for (int idx = top; idx < n; ++idx) {
            int j = xi[idx];
            std::complex<double> pivot = current_pivots[j];
            
            if (std::abs(pivot) < 1e-12) {
                pivot = (pivot.real() < 0.0) ? std::complex<double>(-1e-12, 0.0) : std::complex<double>(1e-12, 0.0);
            }

            std::complex<double> multiplier = x[j] / pivot;
            
            // FIX 2: Decouple multiplier tracking from dense array working space collisions
            l_rows_cols[i].push_back(j);
            l_rows_vals[i].push_back(multiplier);

            for (size_t k_idx = 0; k_idx < u_rows_cols[j].size(); ++k_idx) {
                int col = u_rows_cols[j][k_idx];
                if (col > j) { 
                    x[col] -= multiplier * u_rows_vals[j][k_idx];
                }
            }
        }

        if (std::abs(x[i]) < 1e-12) {
            x[i] = (x[i].real() < 0.0) ? std::complex<double>(-1e-12, 0.0) : std::complex<double>(1e-12, 0.0);
        }
        
        current_pivots[i] = x[i];

        // Store row i's cleared upper-triangular components instantly
        for (int col = i; col < n; ++col) {
            if (std::abs(x[col]) > 1e-15 || col == i) {
                u_rows_cols[i].push_back(col);
                u_rows_vals[i].push_back(x[col]);
            }
            x[col] = std::complex<double>(0.0, 0.0);
        }
        
        // Clean lower-triangular tracking values safely out of the array layout
        for (int idx = top; idx < n; ++idx) {
            x[xi[idx]] = std::complex<double>(0.0, 0.0);
        }
    }

    // Assemble final CSR data matrices for L and U
    L.numRows = n; L.rowPointers.assign(n + 1, 0); L.colIndices.clear(); L.values.clear();
    U.numRows = n; U.rowPointers.assign(n + 1, 0); U.colIndices.clear(); U.values.clear();

    for (int i = 0; i < n; ++i) {
        // FIX 3: Sort L's complex row pairs by column layout index to stay CSR compliant
        std::vector<std::pair<int, std::complex<double>>> l_row_elements;
        for (size_t j = 0; j < l_rows_cols[i].size(); ++j) {
            l_row_elements.push_back({l_rows_cols[i][j], l_rows_vals[i][j]});
        }
        std::sort(l_row_elements.begin(), l_row_elements.end(), 
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        for (const auto& elem : l_row_elements) {
            L.colIndices.push_back(elem.first);
            L.values.push_back(elem.second);
        }
        L.colIndices.push_back(i);
        L.values.push_back(std::complex<double>(1.0, 0.0)); // Unit lower triangular diagonal
        L.rowPointers[i + 1] = static_cast<int>(L.colIndices.size());

        for (size_t j = 0; j < u_rows_cols[i].size(); ++j) {
            U.colIndices.push_back(u_rows_cols[i][j]);
            U.values.push_back(u_rows_vals[i][j]);
        }
        U.rowPointers[i + 1] = static_cast<int>(U.colIndices.size());
    }

    return true;
}
//good and robust
std::vector<std::complex<double>> CircuitMath::luSolveComplex(const CsrMatrixComplex& L, const CsrMatrixComplex& U, const std::vector<int>& P, const std::vector<std::complex<double>>& B) {
    int n = L.numRows;
    // Fix: Explicitly invoke type-safe default constructor signatures
    std::vector<std::complex<double>> x(n, std::complex<double>(0.0, 0.0)), y(n, std::complex<double>(0.0, 0.0));

    // 1. Forward Substitution: L * y = P * B
    for (int i = 0; i < n; ++i) {
        std::complex<double> sum = 0.0;
        int start = L.rowPointers[i];
        int end = L.rowPointers[i + 1]; // Fix: Process total row width securely
        
        for (int j = start; j < end; ++j) {
            int col = L.colIndices[j];
            if (col == i) {
                continue; // Fix: Filter out unit diagonal values safely
            }
            sum += L.values[j] * y[col];
        }
        y[i] = B[P[i]] - sum;
    }

    // 2. Backward Substitution: U * x = y
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = 0.0;
        int start = U.rowPointers[i];
        int end = U.rowPointers[i + 1];
        std::complex<double> diagVal = 1.0;
        
        for (int j = start; j < end; ++j) {
            if (U.colIndices[j] == i) {
                diagVal = U.values[j];
            } else {
                sum += U.values[j] * x[U.colIndices[j]];
            }
        }
        
        if (std::abs(diagVal) < 1e-12) {
            diagVal = (diagVal.real() < 0.0) ? std::complex<double>(-1e-12, 0.0) : std::complex<double>(1e-12, 0.0);
        }
        x[i] = (y[i] - sum) / diagVal;
    }

    return x;
}

// ============================================================================
// ADVANCED DAMPING NON-LINEAR PRIMITIVE LOGIC
// ============================================================================
//good and robust
double CircuitMath::pnjlim(double vnew, double vold, double vt, double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > (2.0 * vt)) {
        if (vold > vcrit) {
            double arg = 1.0 + (vnew - vold) / vt;
            // Fix: Add explicit safety guard against zero/negative logs during extreme updates
            if (arg > 1e-3 && vnew > vold) {
                vnew = vold + vt * std::log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            double arg = (vnew - vcrit) / vt;
            if (arg > -0.99) {
                vnew = vcrit + vt * std::log(1.0 + arg);
            } else {
                vnew = vcrit;
            }
        }
    }
    return vnew;
}
