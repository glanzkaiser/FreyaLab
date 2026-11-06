#include "freyalab/math.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

void MathSolver::accumulateCsr(CsrMatrix& A, int row, int col, double val) 
{
	int start = A.rowPointers[row];
	int end = A.rowPointers[row + 1];

	for (int i = start; i < end; ++i) 
	{
		if (A.colIndices[i] == col) 
		{
			A.values[i] += val;
 			return;
		}
	}
	if (std::abs(val) < 1e-20) return;

	throw std::runtime_error("Structural Failure: Attempted to stamp matrix coordinate (" + 
                             std::to_string(row) + ", " + std::to_string(col) + 
                             ") which was not registered during structural finalization!");
}

void MathSolver::finalizeStructure(CsrMatrix& A, const std::vector<std::pair<int, int>>& nonZeroPatterns) 
{
	int n = A.numRows; 
	A.clear();
	A.numRows = n; 
	A.rowPointers.assign(n + 1, 0);

	std::vector<std::vector<int>> rowCols(n);
	for (const auto& p : nonZeroPatterns) 
	{
		if (p.first >= 0 && p.first < n && p.second >= 0 && p.second < n) 
		{
			rowCols[p.first].push_back(p.second);
		}
	}

	int currentNnz = 0;
	for (int i = 0; i < n; ++i) 
	{
		A.rowPointers[i] = currentNnz;
		std::sort(rowCols[i].begin(), rowCols[i].end());
		rowCols[i].erase(std::unique(rowCols[i].begin(), rowCols[i].end()), rowCols[i].end());

		for (int col : rowCols[i]) 
		{
			A.colIndices.push_back(col);
			A.values.push_back(0.0);
			currentNnz++;
		}
	}
	A.rowPointers[n] = currentNnz;
}

bool MathSolver::sparseLuFactorize(const CsrMatrix& A, CsrMatrix& L, CsrMatrix& U, std::vector<int>& pivot) 
{
	int n = A.numRows;
    
	// Contiguous flat memory allocation replaces nested layout loops
	std::vector<double> denseA(n * n, 0.0);
	for (int r = 0; r < n; ++r) 
	{
		int start = A.rowPointers[r];
		int end = A.rowPointers[r + 1];
		for (int i = start; i < end; ++i) 
		{
			denseA[r * n + A.colIndices[i]] = A.values[i];
		}
	}

	for (int i = 0; i < n; ++i) pivot[i] = i;

	for (int i = 0; i < n; ++i) 
	{
		double maxVal = std::abs(denseA[i * n + i]);
		int pivotRow = i;
		for (int r = i + 1; r < n; ++r) 
		{
			double currentVal = std::abs(denseA[r * n + i]);
			if (currentVal > maxVal) 
			{
				maxVal = currentVal;
				pivotRow = r;
			}
		}

		if (maxVal < 1e-20) return false;

		if (pivotRow != i) 
		{
			std::swap_ranges(denseA.begin() + (i * n), denseA.begin() + ((i + 1) * n), denseA.begin() + (pivotRow * n));
			std::swap(pivot[i], pivot[pivotRow]);
		}

		for (int r = i + 1; r < n; ++r) 
		{
			double factor = denseA[r * n + i] / denseA[i * n + i];
			denseA[r * n + i] = factor; 
			for (int c = i + 1; c < n; ++c) 
			{
				denseA[r * n + c] -= factor * denseA[i * n + c];
			}
		}
	}

	L.clear(); L.numRows = n; L.rowPointers.assign(n + 1, 0);
	U.clear(); U.numRows = n; U.rowPointers.assign(n + 1, 0);

	for (int r = 0; r < n; ++r) 
	{
		for (int c = 0; c < n; ++c) 
		{
			double val = denseA[r * n + c];
			if (c < r) 
			{
				if (std::abs(val) > 1e-20) 
				{
					L.colIndices.push_back(c);
					L.values.push_back(val);
				}
			} 
			else if (c == r) 
			{
				L.colIndices.push_back(c);
				L.values.push_back(1.0); 
				U.colIndices.push_back(c);
				U.values.push_back(val);
			} 
			else 
			{
				if (std::abs(val) > 1e-20) 
				{
					U.colIndices.push_back(c);
					U.values.push_back(val);
				}
			}
		}
		L.rowPointers[r + 1] = static_cast<int>(L.values.size());
		U.rowPointers[r + 1] = static_cast<int>(U.values.size());
	}
	return true;
}

void MathSolver::sparseLuSolve(const CsrMatrix& L, const CsrMatrix& U, const std::vector<int>& pivot, const std::vector<double>& B, std::vector<double>& x) 
{
	int n = L.numRows;
	x.assign(n, 0.0);
	std::vector<double> y(n, 0.0);

	for (int i = 0; i < n; ++i) 
	{
		double sum = 0.0;
		int start = L.rowPointers[i];
		int end = L.rowPointers[i + 1];
		for (int idx = start; idx < end; ++idx) 
		{
			int col = L.colIndices[idx];
			if (col < i) sum += L.values[idx] * y[col];
		}
		y[i] = B[pivot[i]] - sum; 
	}

	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
		double diagVal = 1.0;
		bool diagFound = false;
		int start = U.rowPointers[i];
		int end = U.rowPointers[i + 1];
        
		for (int idx = start; idx < end; ++idx) 
		{
			int col = U.colIndices[idx];
			if (col == i) 
			{
				diagVal = U.values[idx];
				diagFound = true;
			} 
			else if (col > i) 
			{
				sum += U.values[idx] * x[col];
			}
		}
        
		if (!diagFound || std::abs(diagVal) < 1e-20) 
		{
			throw std::runtime_error("Mathematical singularity encountered during backward substitution.");
		}
		x[i] = (y[i] - sum) / diagVal;
	}
}

double MathSolver::pnjlim(double vnew, double vold, double vt, double vcrit) 
{
	if (vnew > vcrit && std::abs(vnew - vold) > (2.0 * vt)) 
	{
		if (vold > 0.0) 
		{
			double arg = 1.0 + (vnew - vold) / vt;
			if (arg > 0.0) 
			{
				vnew = vold + vt * std::log(arg);
			} 
			else 
			{
				vnew = vcrit;
			}
		} 
		else 
		{
			vnew = vt * std::log(vnew / vt);
		}
	} 
	else if (vnew < 0.0) 
	{
		// Standard SPICE PN-junction limiting for negative steps
		double arg = 1.0 + (vnew - vold) / vt;
		if (arg > 0.0) 
		{
			vnew = vold + vt * std::log(arg);
		} 
		else 
		{
			vnew = vold - vt * 2.0;
		}
	}
	return vnew;
}


/*double MathSolver::pnjlim(double vnew, double vold, double vt, double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > (2.0 * vt)) {
        if (vold > 0.0) {
            double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) {
                vnew = vold + vt * std::log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt * std::log(vnew / vt);
        }
    } else if (vnew < 0.0) {
        if (vold > 0.0) {
            // Smooth slope attenuation eliminates non-linear step-cliff drop issues
            double arg = vnew / vt;
            if (arg > -50.0) {
                vnew = vt * std::log(std::exp(arg) + 1e-12);
            } else {
                vnew = vold - vt * 2.0; 
            }
        } else {
            double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) vnew = vold + vt * std::log(arg);
        }
    }
    return vnew;
}*/
