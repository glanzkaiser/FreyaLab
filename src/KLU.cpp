// KLU.cpp
#include "freyalab/KLU.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stack>
#include <utility>

namespace {

// ============================================================================
// 1. MAXIMUM TRANSVERSAL MATCHING (AUGMENTING PATHS) FOR ZERO-FREE DIAGONAL
// ============================================================================
bool findAugmentingPath(int row,
                        const CsrMatrix& A,
                        std::vector<int>& colMatch,
                        std::vector<int>& visitedCol,
                        int epoch) {
    int start = A.rowPointers[row];
    int end   = A.rowPointers[row + 1];

    // First pass: look for an unmatched column immediately (cheap greedy match)
    for (int idx = start; idx < end; ++idx) {
        int col = A.colIndices[idx];
        if (colMatch[col] == -1) {
            colMatch[col] = row;
            return true;
        }
    }

    // Second pass: DFS alternating path
    for (int idx = start; idx < end; ++idx) {
        int col = A.colIndices[idx];
        if (visitedCol[col] == epoch) continue;
        visitedCol[col] = epoch;

        int prevRow = colMatch[col];
        if (prevRow != -1 && findAugmentingPath(prevRow, A, colMatch, visitedCol, epoch)) {
            colMatch[col] = row;
            return true;
        }
    }
    return false;
}

std::vector<int> computeMaxTransversalRowPerm(const CsrMatrix& A) {
    int n = A.numRows;
    std::vector<int> colMatch(n, -1);   // col -> matched original row
    std::vector<int> visitedCol(n, 0);
    int epoch = 1;

    for (int r = 0; r < n; ++r) {
        findAugmentingPath(r, A, colMatch, visitedCol, epoch++);
    }

    // Build rowPerm[newRow] = origRow so that newRow == col has a nonzero diagonal
    std::vector<int> rowPerm(n, -1);
    std::vector<bool> rowUsed(n, false);
    for (int c = 0; c < n; ++c) {
        if (colMatch[c] != -1) {
            rowPerm[c] = colMatch[c];
            rowUsed[colMatch[c]] = true;
        }
    }

    // Fill any structurally singular fallback rows sequentially
    int nextUnused = 0;
    for (int c = 0; c < n; ++c) {
        if (rowPerm[c] == -1) {
            while (nextUnused < n && rowUsed[nextUnused]) ++nextUnused;
            if (nextUnused < n) {
                rowPerm[c] = nextUnused;
                rowUsed[nextUnused] = true;
            } else {
                rowPerm[c] = c;
            }
        }
    }
    return rowPerm;
}

// ============================================================================
// 2. TARJAN STRONGLY CONNECTED COMPONENTS (BTF DECOMPOSITION)
// ============================================================================
void tarjanDfs(int u,
               const std::vector<std::vector<int>>& adj,
               int& timer,
               std::vector<int>& disc,
               std::vector<int>& low,
               std::vector<bool>& onStack,
               std::stack<int>& st,
               std::vector<std::vector<int>>& sccs) {
    disc[u] = low[u] = timer++;
    st.push(u);
    onStack[u] = true;

    for (int v : adj[u]) {
        if (disc[v] == -1) {
            tarjanDfs(v, adj, timer, disc, low, onStack, st, sccs);
            low[u] = std::min(low[u], low[v]);
        } else if (onStack[v]) {
            low[u] = std::min(low[u], disc[v]);
        }
    }

    if (low[u] == disc[u]) {
        std::vector<int> component;
        while (true) {
            int w = st.top();
            st.pop();
            onStack[w] = false;
            component.push_back(w);
            if (w == u) break;
        }
        sccs.push_back(std::move(component));
    }
}

// ============================================================================
// 3. APPROXIMATE MINIMUM DEGREE (AMD) FILL-REDUCING ORDERING PER BLOCK
// ============================================================================
std::vector<int> computeBlockAmdOrdering(int blockSize, const std::vector<std::vector<int>>& localAdj) {
    std::vector<int> order;
    order.reserve(blockSize);
    if (blockSize <= 2) {
        for (int i = 0; i < blockSize; ++i) order.push_back(i);
        return order;
    }

    // Build symmetrized undirected quotient graph G(|B| + |B|^T)
    std::vector<std::vector<int>> symAdj(blockSize);
    for (int u = 0; u < blockSize; ++u) {
        for (int v : localAdj[u]) {
            if (v != u && v >= 0 && v < blockSize) {
                symAdj[u].push_back(v);
                symAdj[v].push_back(u);
            }
        }
    }
    for (int u = 0; u < blockSize; ++u) {
        std::sort(symAdj[u].begin(), symAdj[u].end());
        symAdj[u].erase(std::unique(symAdj[u].begin(), symAdj[u].end()), symAdj[u].end());
    }

    std::vector<bool> eliminated(blockSize, false);
    std::vector<int> degree(blockSize, 0);
    std::vector<int> marker(blockSize, -1);

    // Bucket queue by current approximate degree for O(N) minimum-degree selection
    std::vector<std::vector<int>> buckets(blockSize);
    int minDeg = blockSize;
    for (int u = 0; u < blockSize; ++u) {
        degree[u] = static_cast<int>(symAdj[u].size());
        buckets[degree[u]].push_back(u);
        minDeg = std::min(minDeg, degree[u]);
    }

    for (int step = 0; step < blockSize; ++step) {
        int pivot = -1;
        while (minDeg < blockSize) {
            while (!buckets[minDeg].empty()) {
                int cand = buckets[minDeg].back();
                buckets[minDeg].pop_back();
                if (!eliminated[cand] && degree[cand] == minDeg) {
                    pivot = cand;
                    break;
                }
            }
            if (pivot != -1) break;
            ++minDeg;
        }

        if (pivot == -1) {
            for (int u = 0; u < blockSize; ++u) {
                if (!eliminated[u]) { pivot = u; break; }
            }
        }

        eliminated[pivot] = true;
        order.push_back(pivot);

        // Collect active neighbors of pivot to form fill-in clique
        std::vector<int> neighbors;
        for (int nb : symAdj[pivot]) {
            if (!eliminated[nb]) {
                neighbors.push_back(nb);
            }
        }

        // Update neighbor adjacency & approximate external degrees
        for (size_t i = 0; i < neighbors.size(); ++i) {
            int u = neighbors[i];
            int epochTag = step * blockSize + static_cast<int>(i);

            std::vector<int> newNeighbors;
            newNeighbors.reserve(symAdj[u].size() + neighbors.size());

            for (int w : symAdj[u]) {
                if (!eliminated[w] && marker[w] != epochTag) {
                    marker[w] = epochTag;
                    newNeighbors.push_back(w);
                }
            }
            for (int w : neighbors) {
                if (w != u && !eliminated[w] && marker[w] != epochTag) {
                    marker[w] = epochTag;
                    newNeighbors.push_back(w);
                }
            }

            symAdj[u] = std::move(newNeighbors);
            int newDeg = std::min(blockSize - 1, static_cast<int>(symAdj[u].size()));
            degree[u] = newDeg;
            buckets[newDeg].push_back(u);
            if (newDeg < minDeg) minDeg = newDeg;
        }
    }

    return order;
}

// Sparse topological DFS over previously factored L columns (column-oriented Gilbert-Peierls)
void blockReachDfsCol(int j,
                      const std::vector<std::vector<int>>& lRows,
                      std::vector<int>& visitedEpoch,
                      int epoch,
                      std::vector<int>& xi,
                      int& top) {
    visitedEpoch[j] = epoch;
    for (int row : lRows[j]) {
        if (visitedEpoch[row] != epoch) {
            blockReachDfsCol(row, lRows, visitedEpoch, epoch, xi, top);
        }
    }
    xi[--top] = j;
}


} // namespace

// ============================================================================
// KLU SYMBOLIC ANALYSIS (RUNS ONCE PER CIRCUIT TOPOLOGY)
// ============================================================================
KluSymbolicContext KluSolver::analyzeSymbolic(const CsrMatrix& pattern) {
    KluSymbolicContext sym;
    int n = pattern.numRows;
    sym.n = n;
    if (n <= 0) {
        sym.isAnalyzed = true;
        return sym;
    }

    // Step 1: Maximum Transversal Matching -> ensures nonzero diagonal
    std::vector<int> matchedRowPerm = computeMaxTransversalRowPerm(pattern);

    // Build directed graph on matched vertices: edge (u -> v) if A[matchedRowPerm[u], v] != 0
    std::vector<std::vector<int>> directedAdj(n);
    for (int u = 0; u < n; ++u) {
        int origRow = matchedRowPerm[u];
        int start = pattern.rowPointers[origRow];
        int end   = pattern.rowPointers[origRow + 1];
        for (int idx = start; idx < end; ++idx) {
            int v = pattern.colIndices[idx];
            if (v != u) {
                directedAdj[u].push_back(v);
            }
        }
    }

    // Step 2: Tarjan's SCC produces reverse topological order of BTF blocks
    // Reversing the SCC list puts A into upper block-triangular form!
    int timer = 0;
    std::vector<int> disc(n, -1), low(n, -1);
    std::vector<bool> onStack(n, false);
    std::stack<int> st;
    std::vector<std::vector<int>> sccs;
    for (int i = 0; i < n; ++i) {
        if (disc[i] == -1) {
            tarjanDfs(i, directedAdj, timer, disc, low, onStack, st, sccs);
        }
    }
    std::reverse(sccs.begin(), sccs.end());

    sym.rowPerm.resize(n);
    sym.invRowPerm.resize(n);
    sym.colPerm.resize(n);
    sym.invColPerm.resize(n);

    int cursor = 0;
    sym.maxBlockSize = 1;

    for (const auto& scc : sccs) {
        int bSize = static_cast<int>(scc.size());
        sym.maxBlockSize = std::max(sym.maxBlockSize, bSize);

        // Map local block nodes to compute fill-reducing AMD ordering
        std::vector<int> globalToLocal(n, -1);
        for (int i = 0; i < bSize; ++i) {
            globalToLocal[scc[i]] = i;
        }

        std::vector<std::vector<int>> localAdj(bSize);
        for (int i = 0; i < bSize; ++i) {
            int u = scc[i];
            for (int v : directedAdj[u]) {
                int lv = globalToLocal[v];
                if (lv != -1) {
                    localAdj[i].push_back(lv);
                }
            }
        }

        // Step 3: AMD ordering inside this BTF diagonal block
        std::vector<int> amdOrder = computeBlockAmdOrdering(bSize, localAdj);

        KluBlockFactor blk;
        blk.startIdx = cursor;
        blk.size = bSize;
        blk.amdOrder = amdOrder;
        blk.amdInvOrder.resize(bSize);
        blk.localPivotRow.resize(bSize);
        blk.diagU.assign(bSize, 1.0);
        blk.diagU_c.assign(bSize, std::complex<double>(1.0, 0.0));

        for (int newLocal = 0; newLocal < bSize; ++newLocal) {
            int oldLocal = amdOrder[newLocal];
            blk.amdInvOrder[oldLocal] = newLocal;

            int matchedVert = scc[oldLocal];
            int origRow = matchedRowPerm[matchedVert];
            int origCol = matchedVert;

            sym.rowPerm[cursor] = origRow;
            sym.invRowPerm[origRow] = cursor;
            sym.colPerm[cursor] = origCol;
            sym.invColPerm[origCol] = cursor;
            ++cursor;
        }

        sym.blocks.push_back(std::move(blk));
    }

    // Step 4: Pre-build Off-Diagonal Upper BTF CSR structures (real & complex)
    std::vector<int> vertToBlock(n, 0);
    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        const auto& blk = sym.blocks[bIdx];
        for (int k = 0; k < blk.size; ++k) {
            vertToBlock[blk.startIdx + k] = bIdx;
        }
    }

    sym.offDiagUpper.numRows = n;
    sym.offDiagUpper.rowPointers.assign(n + 1, 0);
    sym.offDiagUpper.colIndices.clear();
    sym.offDiagUpper.values.clear();

    for (int btfRow = 0; btfRow < n; ++btfRow) {
        int origRow = sym.rowPerm[btfRow];
        int myBlock = vertToBlock[btfRow];

        std::vector<int> upperCols;
        int start = pattern.rowPointers[origRow];
        int end   = pattern.rowPointers[origRow + 1];
        for (int idx = start; idx < end; ++idx) {
            int btfCol = sym.invColPerm[pattern.colIndices[idx]];
            if (vertToBlock[btfCol] > myBlock) {
                upperCols.push_back(btfCol);
            }
        }
        std::sort(upperCols.begin(), upperCols.end());
        for (int c : upperCols) {
            sym.offDiagUpper.colIndices.push_back(c);
            sym.offDiagUpper.values.push_back(0.0);
        }
        sym.offDiagUpper.rowPointers[btfRow + 1] = static_cast<int>(sym.offDiagUpper.colIndices.size());
    }
    sym.offDiagUpper_c.initializeFromPattern(sym.offDiagUpper);

    // Pre-allocate zero-allocation real & complex numeric scratch buffers
    sym.workX.assign(sym.maxBlockSize, 0.0);
    sym.visitedEpoch.assign(sym.maxBlockSize, 0);
    sym.touchedEpoch.assign(sym.maxBlockSize, 0);
    sym.touchedCols.reserve(sym.maxBlockSize);
    sym.xi.assign(sym.maxBlockSize, 0);
    sym.btfRhs.assign(n, 0.0);
    sym.btfSol.assign(n, 0.0);

    sym.workX_c.assign(sym.maxBlockSize, std::complex<double>(0.0, 0.0));
    sym.btfRhs_c.assign(n, std::complex<double>(0.0, 0.0));
    sym.btfSol_c.assign(n, std::complex<double>(0.0, 0.0));

    sym.currentEpoch = 1;
    sym.isAnalyzed = true;

    return sym;
}

// ============================================================================
// KLU FAST REAL NUMERIC FACTORIZATION
// ============================================================================
bool KluSolver::factorizeNumeric(const CsrMatrix& A, KluSymbolicContext& sym) {
    if (!sym.isAnalyzed || sym.n != A.numRows) {
        sym = analyzeSymbolic(A);
    }
    int n = sym.n;
    if (n <= 0) return true;

    std::fill(sym.offDiagUpper.values.begin(), sym.offDiagUpper.values.end(), 0.0);
    std::vector<int> vertToBlock(n, 0);
    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        const auto& blk = sym.blocks[bIdx];
        for (int k = 0; k < blk.size; ++k) {
            vertToBlock[blk.startIdx + k] = bIdx;
        }
    }

    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        auto& blk = sym.blocks[bIdx];
        int bStart = blk.startIdx;
        int bSize  = blk.size;

        // Fast O(1) path for 1x1 singleton BTF blocks
        if (bSize == 1) {
            int origRow = sym.rowPerm[bStart];
            int targetOrigCol = sym.colPerm[bStart];
            double diagVal = 0.0;

            int rStart = A.rowPointers[origRow];
            int rEnd   = A.rowPointers[origRow + 1];
            int offPtr = sym.offDiagUpper.rowPointers[bStart];
            int offEnd = sym.offDiagUpper.rowPointers[bStart + 1];

            for (int idx = rStart; idx < rEnd; ++idx) {
                int origCol = A.colIndices[idx];
                double val = A.values[idx];
                if (origCol == targetOrigCol) {
                    diagVal = val;
                } else {
                    int btfCol = sym.invColPerm[origCol];
                    if (vertToBlock[btfCol] > bIdx) {
                        auto it = std::lower_bound(sym.offDiagUpper.colIndices.begin() + offPtr,
                                                   sym.offDiagUpper.colIndices.begin() + offEnd,
                                                   btfCol);
                        if (it != sym.offDiagUpper.colIndices.begin() + offEnd && *it == btfCol) {
                            sym.offDiagUpper.values[it - sym.offDiagUpper.colIndices.begin()] = val;
                        }
                    }
                }
            }

            if (std::abs(diagVal) < 1e-14) {
                diagVal = (diagVal < 0.0) ? -1e-14 : 1e-14;
            }
            blk.diagU[0] = diagVal;
            blk.localPivotRow[0] = 0;
            continue;
        }

        // General multi-node BTF block sparse column-oriented left-looking LU with Partial Pivoting
        constexpr double PIVOT_TOL = 0.1;

        std::vector<std::vector<std::pair<int, double>>> colEntries(bSize);
        for (int localRow = 0; localRow < bSize; ++localRow) {
            int btfRow = bStart + localRow;
            int origRow = sym.rowPerm[btfRow];

            int rStart   = A.rowPointers[origRow];
            int rEnd     = A.rowPointers[origRow + 1];
            int offStart = sym.offDiagUpper.rowPointers[btfRow];
            int offEnd   = sym.offDiagUpper.rowPointers[btfRow + 1];

            for (int idx = rStart; idx < rEnd; ++idx) {
                int btfCol = sym.invColPerm[A.colIndices[idx]];
                double val = A.values[idx];

                if (btfCol >= bStart && btfCol < bStart + bSize) {
                    int localCol = btfCol - bStart;
                    colEntries[localCol].push_back({localRow, val});
                } else if (vertToBlock[btfCol] > bIdx) {
                    auto it = std::lower_bound(sym.offDiagUpper.colIndices.begin() + offStart,
                                               sym.offDiagUpper.colIndices.begin() + offEnd,
                                               btfCol);
                    if (it != sym.offDiagUpper.colIndices.begin() + offEnd && *it == btfCol) {
                        sym.offDiagUpper.values[it - sym.offDiagUpper.colIndices.begin()] = val;
                    }
                }
            }
        }

        std::vector<std::vector<int>> lRows(bSize), uRows(bSize);
        std::vector<std::vector<double>> lVals(bSize), uVals(bSize);
        std::vector<int> rowToStep(bSize, -1);
        blk.localPivotRow.resize(bSize);

        for (int k = 0; k < bSize; ++k) {
            int epoch = ++sym.currentEpoch;
            int top = bSize;
            sym.touchedCols.clear();

            auto markAndAccumulate = [&](int row, double val) {
                if (sym.touchedEpoch[row] != epoch) {
                    sym.touchedEpoch[row] = epoch;
                    sym.workX[row] = val;
                    sym.touchedCols.push_back(row);
                } else {
                    sym.workX[row] += val;
                }
            };

            // 1. Scatter column k entries & find topological reach across already-pivoted rows
            for (const auto& entry : colEntries[k]) {
                int r = entry.first;
                markAndAccumulate(r, entry.second);
                int step = rowToStep[r];
                if (step != -1 && sym.visitedEpoch[step] != epoch) {
                    blockReachDfsCol(step, lRows, sym.visitedEpoch, epoch, sym.xi, top);
                }
            }

            // 2. Triangular solve L(:, 0..k-1) * x = a_k in topological order
            for (int idx = top; idx < bSize; ++idx) {
                int j = sym.xi[idx];
                int pivRow = blk.localPivotRow[j];
                double xj = (sym.touchedEpoch[pivRow] == epoch) ? sym.workX[pivRow] : 0.0;

                if (std::abs(xj) > 1e-18) {
                    for (size_t p = 0; p < lRows[j].size(); ++p) {
                        markAndAccumulate(lRows[j][p], -xj * lVals[j][p]);
                    }
                }
            }

            // 3. Threshold Partial Pivoting among unpivoted rows (rowToStep[r] == -1)
            int bestRow = -1;
            double maxAbs = 0.0;
            double diagCandAbs = 0.0;
            int diagCandRow = -1;

            for (int r : sym.touchedCols) {
                if (rowToStep[r] == -1) {
                    double absVal = std::abs(sym.workX[r]);
                    if (absVal > maxAbs) {
                        maxAbs = absVal;
                        bestRow = r;
                    }
                    if (r == k) {
                        diagCandAbs = absVal;
                        diagCandRow = r;
                    }
                }
            }

            if (diagCandRow != -1 && diagCandAbs >= PIVOT_TOL * maxAbs && diagCandAbs > 1e-14) {
                bestRow = diagCandRow;
            }
            if (bestRow == -1) {
                for (int r = 0; r < bSize; ++r) {
                    if (rowToStep[r] == -1) { bestRow = r; break; }
                }
            }

            rowToStep[bestRow] = k;
            blk.localPivotRow[k] = bestRow;

            double pivotVal = (sym.touchedEpoch[bestRow] == epoch) ? sym.workX[bestRow] : 0.0;
            if (std::abs(pivotVal) < 1e-14) {
                pivotVal = (pivotVal < 0.0) ? -1e-14 : 1e-14;
            }
            blk.diagU[k] = pivotVal;

            // 4. Separate U(:, k) (pivoted rows step < k) and L(:, k) (unpivoted rows / pivotVal)
            for (int r : sym.touchedCols) {
                double val = sym.workX[r];
                if (r != bestRow && std::abs(val) > 1e-18) {
                    int step = rowToStep[r];
                    if (step != -1 && step < k) {
                        uRows[k].push_back(step);
                        uVals[k].push_back(val);
                    } else if (step == -1) {
                        lRows[k].push_back(r);
                        lVals[k].push_back(val / pivotVal);
                    }
                }
                sym.workX[r] = 0.0;
            }
        }

        // 5. Pack L and U into row-major CSR for fast block substitution in solveNumeric
        std::vector<std::vector<std::pair<int, double>>> rowL(bSize), rowU(bSize);
        for (int col = 0; col < bSize; ++col) {
            for (size_t p = 0; p < lRows[col].size(); ++p) {
                int pivStepRow = rowToStep[lRows[col][p]];
                rowL[pivStepRow].push_back({col, lVals[col][p]});
            }
            for (size_t p = 0; p < uRows[col].size(); ++p) {
                int pivStepRow = uRows[col][p];
                rowU[pivStepRow].push_back({col, uVals[col][p]});
            }
        }

        blk.L.numRows = bSize;
        blk.L.rowPointers.assign(bSize + 1, 0);
        blk.L.colIndices.clear();
        blk.L.values.clear();

        blk.U.numRows = bSize;
        blk.U.rowPointers.assign(bSize + 1, 0);
        blk.U.colIndices.clear();
        blk.U.values.clear();

        for (int i = 0; i < bSize; ++i) {
            std::sort(rowL[i].begin(), rowL[i].end());
            for (const auto& kv : rowL[i]) {
                blk.L.colIndices.push_back(kv.first);
                blk.L.values.push_back(kv.second);
            }
            blk.L.rowPointers[i + 1] = static_cast<int>(blk.L.colIndices.size());

            std::sort(rowU[i].begin(), rowU[i].end());
            for (const auto& kv : rowU[i]) {
                blk.U.colIndices.push_back(kv.first);
                blk.U.values.push_back(kv.second);
            }
            blk.U.rowPointers[i + 1] = static_cast<int>(blk.U.colIndices.size());
        }
    }

    return true;
}

// ============================================================================
// KLU REAL BLOCK SUBSTITUTION SOLVE
// ============================================================================
std::vector<double> KluSolver::solveNumeric(KluSymbolicContext& sym, const std::vector<double>& b) {
    int n = sym.n;
    std::vector<double> x(n, 0.0);
    if (n <= 0 || static_cast<int>(b.size()) < n) return x;

    for (int btfRow = 0; btfRow < n; ++btfRow) {
        sym.btfRhs[btfRow] = b[sym.rowPerm[btfRow]];
    }

    for (int bIdx = static_cast<int>(sym.blocks.size()) - 1; bIdx >= 0; --bIdx) {
        const auto& blk = sym.blocks[bIdx];
        int bStart = blk.startIdx;
        int bSize  = blk.size;

        for (int i = 0; i < bSize; ++i) {
            int btfRow = bStart + i;
            double sumUpper = 0.0;
            int uStart = sym.offDiagUpper.rowPointers[btfRow];
            int uEnd   = sym.offDiagUpper.rowPointers[btfRow + 1];
            for (int idx = uStart; idx < uEnd; ++idx) {
                sumUpper += sym.offDiagUpper.values[idx] * sym.btfSol[sym.offDiagUpper.colIndices[idx]];
            }
            sym.btfRhs[btfRow] -= sumUpper;
        }

        if (bSize == 1) {
            sym.btfSol[bStart] = sym.btfRhs[bStart] / blk.diagU[0];
            continue;
        }

        // Forward substitution inside block using intra-block row pivot permutation: L * y = P_blk * rhs_blk
        for (int i = 0; i < bSize; ++i) {
            double sum = 0.0;
            int lStart = blk.L.rowPointers[i];
            int lEnd   = blk.L.rowPointers[i + 1];
            for (int idx = lStart; idx < lEnd; ++idx) {
                sum += blk.L.values[idx] * sym.btfSol[bStart + blk.L.colIndices[idx]];
            }
            int pivRow = blk.localPivotRow[i];
            sym.btfSol[bStart + i] = sym.btfRhs[bStart + pivRow] - sum;
        }


        for (int i = bSize - 1; i >= 0; --i) {
            double sum = 0.0;
            int uStart = blk.U.rowPointers[i];
            int uEnd   = blk.U.rowPointers[i + 1];
            for (int idx = uStart; idx < uEnd; ++idx) {
                sum += blk.U.values[idx] * sym.btfSol[bStart + blk.U.colIndices[idx]];
            }
            sym.btfSol[bStart + i] = (sym.btfSol[bStart + i] - sum) / blk.diagU[i];
        }
    }

    for (int btfCol = 0; btfCol < n; ++btfCol) {
        x[sym.colPerm[btfCol]] = sym.btfSol[btfCol];
    }

    return x;
}

// ============================================================================
// KLU FAST COMPLEX NUMERIC FACTORIZATION (AC PHASOR / BODE SWEEPS)
// ============================================================================
bool KluSolver::factorizeNumericComplex(const CsrMatrixComplex& A, KluSymbolicContext& sym) {
    int n = sym.n;
    if (!sym.isAnalyzed || n != A.numRows || n <= 0) return false;

    // 1. Inisialisasi lazy workspace kompleks jika belum dialokasikan
    if (static_cast<int>(sym.offDiagUpper_c.values.size()) != static_cast<int>(sym.offDiagUpper.colIndices.size())) {
        sym.offDiagUpper_c.initializeFromPattern(sym.offDiagUpper);
        sym.workX_c.assign(sym.maxBlockSize, std::complex<double>(0.0, 0.0));
        sym.btfRhs_c.assign(n, std::complex<double>(0.0, 0.0));
        sym.btfSol_c.assign(n, std::complex<double>(0.0, 0.0));
    }

    std::fill(sym.offDiagUpper_c.values.begin(), sym.offDiagUpper_c.values.end(), std::complex<double>(0.0, 0.0));

    // Pemetaan indeks BTF global ke indeks blok BTF
    std::vector<int> vertToBlock(n, 0);
    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        const auto& blk = sym.blocks[bIdx];
        for (int k = 0; k < blk.size; ++k) {
            vertToBlock[blk.startIdx + k] = bIdx;
        }
    }

    // 2. Ekstraksi seluruh blok BTF dan off-diagonal upper secara column-oriented
    //    Struktur per blok: blockCols[bIdx][amdCol] berisi pasangan (amdRow, val)
    struct EntryC {
        int row;
        std::complex<double> val;
    };
    std::vector<std::vector<std::vector<EntryC>>> blockCols(sym.blocks.size());
    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        blockCols[bIdx].resize(sym.blocks[bIdx].size);
    }

    for (int origRow = 0; origRow < n; ++origRow) {
        int btfRow = sym.invRowPerm[origRow];
        int rBlk   = vertToBlock[btfRow];

        int rStart   = A.rowPointers[origRow];
        int rEnd     = A.rowPointers[origRow + 1];
        int offStart = sym.offDiagUpper_c.rowPointers[btfRow];
        int offEnd   = sym.offDiagUpper_c.rowPointers[btfRow + 1];

        for (int idx = rStart; idx < rEnd; ++idx) {
            int btfCol = sym.invColPerm[A.colIndices[idx]];
            int cBlk   = vertToBlock[btfCol];
            std::complex<double> val = A.values[idx];

            if (rBlk == cBlk) {
                const auto& blk = sym.blocks[rBlk];
                int rawRow = btfRow - blk.startIdx;
                int rawCol = btfCol - blk.startIdx;
                int amdRow = blk.amdInvOrder.empty() ? rawRow : blk.amdInvOrder[rawRow];
                int amdCol = blk.amdInvOrder.empty() ? rawCol : blk.amdInvOrder[rawCol];
                blockCols[rBlk][amdCol].push_back({amdRow, val});
            } else if (cBlk > rBlk) {
                auto it = std::lower_bound(sym.offDiagUpper_c.colIndices.begin() + offStart,
                                           sym.offDiagUpper_c.colIndices.begin() + offEnd,
                                           btfCol);
                if (it != sym.offDiagUpper_c.colIndices.begin() + offEnd && *it == btfCol) {
                    sym.offDiagUpper_c.values[it - sym.offDiagUpper_c.colIndices.begin()] = val;
                }
            }
        }
    }

    // 3. Faktorisasi LU Column-Oriented Gilbert-Peierls per Blok BTF
    for (int bIdx = 0; bIdx < static_cast<int>(sym.blocks.size()); ++bIdx) {
        auto& blk = sym.blocks[bIdx];
        int bSize = blk.size;

        blk.diagU_c.assign(bSize, std::complex<double>(1.0, 0.0));
        blk.localPivotRow.resize(bSize);
        for (int i = 0; i < bSize; ++i) {
            blk.localPivotRow[i] = i;
        }

        // Kasus khusus blok skalar 1x1
        if (bSize == 1) {
            std::complex<double> diagVal(0.0, 0.0);
            for (const auto& entry : blockCols[bIdx][0]) {
                diagVal += entry.val;
            }
            if (std::abs(diagVal) < 1e-14) {
                diagVal = (diagVal.real() < 0.0) ? std::complex<double>(-1e-14, 0.0) : std::complex<double>(1e-14, 0.0);
            }
            blk.diagU_c[0] = diagVal;
            blk.L_c.numRows = 1;
            blk.L_c.rowPointers = {0, 0};
            blk.L_c.colIndices.clear();
            blk.L_c.values.clear();
            blk.U_c.numRows = 1;
            blk.U_c.rowPointers = {0, 0};
            blk.U_c.colIndices.clear();
            blk.U_c.values.clear();
            continue;
        }

        // lRows[j] menyimpan indeks baris (sudah ter-pivot) untuk kolom j dari L
        std::vector<std::vector<int>> lRows(bSize), uRows(bSize);
        std::vector<std::vector<std::complex<double>>> lVals(bSize), uVals(bSize);

        // Permutasi pivot dinamis: rowToStep[origAmdRow] = langkah eliminasi k
        std::vector<int> rowToStep(bSize, -1);

        for (int k = 0; k < bSize; ++k) {
            int epoch = ++sym.currentEpoch;
            if (epoch >= 1000000000) {
                std::fill(sym.visitedEpoch.begin(), sym.visitedEpoch.end(), 0);
                std::fill(sym.touchedEpoch.begin(), sym.touchedEpoch.end(), 0);
                sym.currentEpoch = 1;
                epoch = 1;
            }

            int top = bSize;
            sym.touchedCols.clear();

            auto markAndAccumulate = [&](int row, std::complex<double> val) {
                if (sym.touchedEpoch[row] != epoch) {
                    sym.touchedEpoch[row] = epoch;
                    sym.workX_c[row] = val;
                    sym.touchedCols.push_back(row);
                } else {
                    sym.workX_c[row] += val;
                }
            };

            // Scatter kolom ke-k dari A ke dalam workspace workX_c
            for (const auto& entry : blockCols[bIdx][k]) {
                int r = entry.row;
                markAndAccumulate(r, entry.val);

                // Jika baris r sudah dipivot pada langkah j < k, telusuri reachability kolom L
                int j = rowToStep[r];
                if (j != -1 && sym.visitedEpoch[j] != epoch) {
                    blockReachDfsCol(j, lRows, sym.visitedEpoch, epoch, sym.xi, top);
                }
            }

            // Triangular solve L * u_k = a_k sesuai urutan topologis dari stack xi
            for (int idx = top; idx < bSize; ++idx) {
                int j = sym.xi[idx];
                int pivOrigRow = blk.localPivotRow[j];
                std::complex<double> u_jk = sym.workX_c[pivOrigRow];

                if (std::abs(u_jk) > 1e-18) {
                    uRows[k].push_back(j);
                    uVals[k].push_back(u_jk);

                    for (size_t p = 0; p < lRows[j].size(); ++p) {
                        int targetOrigRow = lRows[j][p];
                        markAndAccumulate(targetOrigRow, -lVals[j][p] * u_jk);
                    }
                }
            }

            // Partial pivoting pada baris yang belum terpilih (rowToStep[r] == -1)
            int bestOrigRow = -1;
            double maxAbs = -1.0;
            for (int r : sym.touchedCols) {
                if (rowToStep[r] == -1) {
                    double mag = std::abs(sym.workX_c[r]);
                    if (mag > maxAbs) {
                        maxAbs = mag;
                        bestOrigRow = r;
                    }
                }
            }

            // Fallback jika semua baris yang tersentuh bernilai 0 (matriks singular struktural)
            if (bestOrigRow == -1) {
                for (int r = 0; r < bSize; ++r) {
                    if (rowToStep[r] == -1) {
                        bestOrigRow = r;
                        break;
                    }
                }
            }

            // Catat pivot untuk langkah k
            blk.localPivotRow[k] = bestOrigRow;
            rowToStep[bestOrigRow] = k;

            std::complex<double> diagVal = (sym.touchedEpoch[bestOrigRow] == epoch)
                                               ? sym.workX_c[bestOrigRow]
                                               : std::complex<double>(0.0, 0.0);

            if (std::abs(diagVal) < 1e-14) {
                diagVal = (diagVal.real() < 0.0) ? std::complex<double>(-1e-14, 0.0)
                                                 : std::complex<double>(1e-14, 0.0);
            }
            blk.diagU_c[k] = diagVal;

            // Normalisasi kolom ke-k dari L (semua baris yang belum dipivot selain bestOrigRow)
            for (int r : sym.touchedCols) {
                if (rowToStep[r] == -1) {
                    std::complex<double> l_rk = sym.workX_c[r] / diagVal;
                    if (std::abs(l_rk) > 1e-18) {
                        lRows[k].push_back(r);
                        lVals[k].push_back(l_rk);
                    }
                }
                sym.workX_c[r] = std::complex<double>(0.0, 0.0);
            }
        }

        // 4. Simpan faktor L_c dan U_c ke format CSR per-kolom
        blk.L_c.numRows = bSize;
        blk.L_c.rowPointers.assign(bSize + 1, 0);
        blk.L_c.colIndices.clear();
        blk.L_c.values.clear();

        blk.U_c.numRows = bSize;
        blk.U_c.rowPointers.assign(bSize + 1, 0);
        blk.U_c.colIndices.clear();
        blk.U_c.values.clear();

        for (int k = 0; k < bSize; ++k) {
            for (size_t p = 0; p < lRows[k].size(); ++p) {
                // Simpan indeks step pivot-nya atau indeks baris lokal sesuai konvensi solveNumericComplex
                blk.L_c.colIndices.push_back(rowToStep[lRows[k][p]]);
                blk.L_c.values.push_back(lVals[k][p]);
            }
            blk.L_c.rowPointers[k + 1] = static_cast<int>(blk.L_c.colIndices.size());

            for (size_t p = 0; p < uRows[k].size(); ++p) {
                blk.U_c.colIndices.push_back(uRows[k][p]);
                blk.U_c.values.push_back(uVals[k][p]);
            }
            blk.U_c.rowPointers[k + 1] = static_cast<int>(blk.U_c.colIndices.size());
        }
    }

    return true;
}


// ============================================================================
// KLU COMPLEX BLOCK SUBSTITUTION SOLVE
// ============================================================================
std::vector<std::complex<double>> KluSolver::solveNumericComplex(KluSymbolicContext& sym,
                                                                   const std::vector<std::complex<double>>& b) {
    int n = sym.n;
    std::vector<std::complex<double>> x(n, std::complex<double>(0.0, 0.0));
    if (n <= 0 || static_cast<int>(b.size()) < n) return x;

    // 1. Permute RHS vector into BTF row order
    for (int btfRow = 0; btfRow < n; ++btfRow) {
        sym.btfRhs_c[btfRow] = b[sym.rowPerm[btfRow]];
    }

    // 2. Solve diagonal BTF blocks in reverse topological order (bottom to top)
    for (int bIdx = static_cast<int>(sym.blocks.size()) - 1; bIdx >= 0; --bIdx) {
        const auto& blk = sym.blocks[bIdx];
        int bStart = blk.startIdx;
        int bSize  = blk.size;

        // Subtract known contributions from already-solved upper blocks: rhs_i -= U_offdiag * x_upper
        for (int i = 0; i < bSize; ++i) {
            int btfRow = bStart + i;
            std::complex<double> sumUpper(0.0, 0.0);
            int uStart = sym.offDiagUpper_c.rowPointers[btfRow];
            int uEnd   = sym.offDiagUpper_c.rowPointers[btfRow + 1];
            for (int idx = uStart; idx < uEnd; ++idx) {
                sumUpper += sym.offDiagUpper_c.values[idx] * sym.btfSol_c[sym.offDiagUpper_c.colIndices[idx]];
            }
            sym.btfRhs_c[btfRow] -= sumUpper;
        }

        // Fast O(1) division for 1x1 singleton BTF blocks
        if (bSize == 1) {
            sym.btfSol_c[bStart] = sym.btfRhs_c[bStart] / blk.diagU_c[0];
            continue;
        }

        // Forward substitution inside block: L_c * y = rhs_blk
        for (int i = 0; i < bSize; ++i) {
            std::complex<double> sum(0.0, 0.0);
            int lStart = blk.L_c.rowPointers[i];
            int lEnd   = blk.L_c.rowPointers[i + 1];
            for (int idx = lStart; idx < lEnd; ++idx) {
                sum += blk.L_c.values[idx] * sym.btfSol_c[bStart + blk.L_c.colIndices[idx]];
            }
            sym.btfSol_c[bStart + i] = sym.btfRhs_c[bStart + i] - sum;
        }

        // Backward substitution inside block: U_c * x_blk = y
        for (int i = bSize - 1; i >= 0; --i) {
            std::complex<double> sum(0.0, 0.0);
            int uStart = blk.U_c.rowPointers[i];
            int uEnd   = blk.U_c.rowPointers[i + 1];
            for (int idx = uStart; idx < uEnd; ++idx) {
                sum += blk.U_c.values[idx] * sym.btfSol_c[bStart + blk.U_c.colIndices[idx]];
            }
            sym.btfSol_c[bStart + i] = (sym.btfSol_c[bStart + i] - sum) / blk.diagU_c[i];
        }
    }

    // 3. Unpermute complex solution back to original MNA column order
    for (int btfCol = 0; btfCol < n; ++btfCol) {
        x[sym.colPerm[btfCol]] = sym.btfSol_c[btfCol];
    }

    return x;
}
