#include "freyalab/BSIM.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>

#if defined(_WIN32)
    #define NOMINMAX
    #include <windows.h>
#else
    #include <dlfcn.h>
#endif

namespace {
    constexpr double EPS_OX_SIO2    = 3.9 * 8.854187817e-12; // F/m
    constexpr double V_THERMAL_300K = 0.025851991;           // kT/q at 300 K (V)
    constexpr double GMIN_FLOOR     = 1e-12;                 // Minimum conductance (S)
    constexpr double TR_GAMMA       = 2.0 - 1.4142135623730951; // TR-BDF2 gamma = 2 - sqrt(2)

    inline void safeSoftplus(double x, double& sp, double& dsp_dx) {
        if (x > 40.0) {
            sp = x;
            dsp_dx = 1.0;
        } else if (x < -40.0) {
            sp = std::exp(x);
            dsp_dx = sp;
        } else {
            double ex = std::exp(x);
            sp = std::log1p(ex);
            dsp_dx = ex / (1.0 + ex);
        }
    }

    std::string toUpperStr(const std::string& s) {
        std::string out = s;
        std::transform(out.begin(), out.end(), out.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        return out;
    }

    // Helper to construct OSDI v0.4 simulation parameters ($simparam)
    struct DefaultOsdiSimParas {
        std::array<char, 5> gminName{'g', 'm', 'i', 'n', '\0'};
        std::array<char*, 2> names{gminName.data(), nullptr};
        std::array<double, 1> vals{GMIN_FLOOR};
        std::array<char*, 1> emptyStr{nullptr};
        OsdiSimParas paras{};

        DefaultOsdiSimParas() {
            paras.names     = names.data();
            paras.vals      = vals.data();
            paras.names_str = emptyStr.data();
            paras.vals_str  = emptyStr.data();
        }
    };

    // Binds OSDI v0.4 resistive & reactive Jacobian double* pointers inside instBuf
    void bindOsdiJacobianSlots(const OsdiDescriptor* desc,
                               uint8_t* instBuf,
                               std::vector<double>& resistSlots,
                               std::vector<double>& reactSlots) {
        const uint32_t nJac = desc->num_jacobian_entries;
        resistSlots.assign(nJac, 0.0);
        reactSlots.assign(nJac, 0.0);

        auto** resistPtrArray = reinterpret_cast<double**>(instBuf + desc->jacobian_ptr_resist_offset);
        for (uint32_t e = 0; e < nJac; ++e) {
            const OsdiJacobianEntry& entry = *(desc->jacobian_entries + e);
            *(resistPtrArray + e) = &resistSlots.at(e);
            if (entry.flags & JACOBIAN_ENTRY_REACT) {
                auto** reactSlotPtr = reinterpret_cast<double**>(instBuf + entry.react_ptr_off);
                *reactSlotPtr = &reactSlots.at(e);
            }
        }
    }
}

// ============================================================================
// BsimCmgOsdiLoader Implementation
// ============================================================================

BsimCmgOsdiLoader::~BsimCmgOsdiLoader() {
    unload();
}

BsimCmgOsdiLoader& BsimCmgOsdiLoader::getGlobalInstance() {
    static BsimCmgOsdiLoader instance;
    if (!instance.attemptedLoad) {
        instance.attemptedLoad = true;
        if (!instance.load("models/bsimcmg.osdi") && !instance.load("bsimcmg.osdi")) {
            std::cout << "[FreyaLab BSIM] 'models/bsimcmg.osdi' not found. "
                      << "Using built-in Phase-1 BSIM-CMG analytical physics engine.\n";
        }
    }
    return instance;
}

bool BsimCmgOsdiLoader::load(const std::string& osdiPath) {
    unload();

#if defined(_WIN32)
    HMODULE mod = LoadLibraryA(osdiPath.c_str());
    if (!mod) return false;
    libHandle = static_cast<void*>(mod);

    auto* numDescPtr = reinterpret_cast<const uint32_t*>(GetProcAddress(mod, "OSDI_NUM_DESCRIPTORS"));
    auto* descArray  = reinterpret_cast<const OsdiDescriptor*>(GetProcAddress(mod, "OSDI_DESCRIPTORS"));
#else
    void* mod = dlopen(osdiPath.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!mod) return false;
    libHandle = mod;

    auto* numDescPtr = reinterpret_cast<const uint32_t*>(dlsym(mod, "OSDI_NUM_DESCRIPTORS"));
    auto* descArray  = reinterpret_cast<const OsdiDescriptor*>(dlsym(mod, "OSDI_DESCRIPTORS"));
#endif

    if (!numDescPtr || !descArray || *numDescPtr == 0) {
        unload();
        return false;
    }

    // Locate the bsimcmg descriptor inside the OSDI descriptor table
    for (uint32_t i = 0; i < *numDescPtr; ++i) {
        const OsdiDescriptor* cand = descArray + i;
        std::string modName = cand->name ? toUpperStr(cand->name) : "";
        if (modName.find("BSIMCMG") != std::string::npos || *numDescPtr == 1) {
            descriptor = cand;
            break;
        }
    }

    if (!descriptor) {
        unload();
        return false;
    }

    std::cout << "[FreyaLab BSIM] Loaded OSDI compiled model '" << descriptor->name
              << "' from '" << osdiPath << "' (" << descriptor->num_params << " parameters).\n";
    return true;
}

void BsimCmgOsdiLoader::unload() {
    if (libHandle) {
#if defined(_WIN32)
        FreeLibrary(static_cast<HMODULE>(libHandle));
#else
        dlclose(libHandle);
#endif
        libHandle = nullptr;
    }
    descriptor = nullptr;
}

// ============================================================================
// BsimCmgFet Constructor & OSDI Setup (OSDI v0.4)
// ============================================================================

BsimCmgFet::BsimCmgFet(const std::string& name,
                       const std::vector<NetNode>& nodes,
                       ComponentType subType,
                       const BsimCmgParameters& modelParams)
    : Component(name, subType, nodes), params(modelParams) {
    if (subType == ComponentType::GAA_FET || subType == ComponentType::RIBBONFET) {
        params.geomod = MultiGateGeometry::GAA_NANOSHEET;
    }

    const int safeNfin = std::max(1, params.NFIN);
    const double safeEot = std::max(params.EOT, 1e-10);
    cox = EPS_OX_SIO2 / safeEot;

    if (params.geomod == MultiGateGeometry::GAA_NANOSHEET) {
        weff = static_cast<double>(safeNfin) * 2.0 * (params.WSHEET + params.TFIN);
    } else {
        weff = static_cast<double>(safeNfin) * (2.0 * params.HFIN + params.TFIN);
    }

    // Attempt binding to compiled BSIM-CMG 112.1.0 OSDI library
    BsimCmgOsdiLoader& loader = BsimCmgOsdiLoader::getGlobalInstance();
    if (loader.isLoaded()) {
        osdiDesc = loader.getDescriptor();
        setupOsdiBuffers();
    }
}

void BsimCmgFet::setOsdiParamByName(const std::string& paramName, double val) {
    if (!osdiDesc) return;
    std::string target = toUpperStr(paramName);

    for (uint32_t i = 0; i < osdiDesc->num_params; ++i) {
        const OsdiParamOpvar& p = *(osdiDesc->param_opvar + i);
        uint32_t kind = p.flags & PARA_KIND_MASK;
        if (kind != PARA_KIND_MODEL && kind != PARA_KIND_INST) {
            continue;
        }

        bool matched = false;
        for (uint32_t a = 0; a < p.num_alias + 1; ++a) {
            const char* aliasPtr = *(p.name + a);
            if (aliasPtr && toUpperStr(aliasPtr) == target) {
                matched = true;
                break;
            }
        }
        if (!matched) continue;

        bool isInstance = (kind == PARA_KIND_INST);
        void* rawPtr = isInstance
            ? osdiDesc->access(instanceData.data(), modelData.data(), i, ACCESS_FLAG_SET | ACCESS_FLAG_INSTANCE)
            : osdiDesc->access(nullptr, modelData.data(), i, ACCESS_FLAG_SET);

        if (!rawPtr) continue;

        uint32_t baseType = p.flags & PARA_TY_MASK;
        if (baseType == PARA_TY_REAL) {
            *reinterpret_cast<double*>(rawPtr) = val;
        } else if (baseType == PARA_TY_INT) {
            *reinterpret_cast<int32_t*>(rawPtr) = static_cast<int32_t>(std::llround(val));
        }
    }
}

void BsimCmgFet::setupOsdiBuffers() {
    if (!osdiDesc) return;

    modelData.assign(osdiDesc->model_size, 0);
    instanceData.assign(osdiDesc->instance_size, 0);

    // 1. Apply default & geometry parameters before setup_model / setup_instance
    setOsdiParamByName("TYPE",    params.isNChannel ? 1.0 : -1.0);
    setOsdiParamByName("GEOMOD",  static_cast<double>(static_cast<int>(params.geomod)));
    setOsdiParamByName("L",       params.L);
    setOsdiParamByName("TFIN",    params.TFIN);
    setOsdiParamByName("HFIN",    params.HFIN);
    setOsdiParamByName("WSHEET",  params.WSHEET);
    setOsdiParamByName("NFIN",    static_cast<double>(params.NFIN));
    setOsdiParamByName("EOT",     params.EOT);
    setOsdiParamByName("VTH0",    params.VTH0);
    setOsdiParamByName("NFACTOR", params.NFACTOR);
    setOsdiParamByName("ETA0",    params.ETA0);
    setOsdiParamByName("U0",      params.U0);
    setOsdiParamByName("VSAT",    params.VSAT);
    setOsdiParamByName("PCLM",    params.PCLM);

    for (const auto& kv : params.rawOverrides) {
        setOsdiParamByName(kv.first, kv.second);
    }

    DefaultOsdiSimParas simEnv;

    // 2. Initialize OSDI v0.4 Model
    OsdiInitInfo modelInit{};
    osdiDesc->setup_model(nullptr, modelData.data(), &simEnv.paras, &modelInit);
    if (modelInit.flags & EVAL_RET_FLAG_FATAL) {
        osdiReady = false;
        return;
    }

    // 3. Configure external terminal mapping (D, G, S, [B])
    uint32_t numTerms = osdiDesc->num_terminals;
    uint32_t connectedTerms = std::min<uint32_t>(numTerms, std::max<uint32_t>(static_cast<uint32_t>(nodes.size()), 4u));
    auto* nodeMappingPtr = reinterpret_cast<uint32_t*>(instanceData.data() + osdiDesc->node_mapping_offset);
    for (uint32_t t = 0; t < numTerms; ++t) {
        *(nodeMappingPtr + t) = (t < nodes.size() || t == 3) ? t : UINT32_MAX;
    }

        // 4. Initialize OSDI v0.4 Instance & Temperature Setup (300.15 K)
    OsdiInitInfo instInit{};
    osdiDesc->setup_instance(nullptr, instanceData.data(), modelData.data(),
                             300.15, connectedTerms, &simEnv.paras, &instInit);
    if (instInit.flags & EVAL_RET_FLAG_FATAL) {
        osdiReady = false;
        return;
    }

    // Map internal node indices & state indices
    uint32_t totalOsdiNodes = osdiDesc->num_nodes;
    for (uint32_t n = numTerms; n < totalOsdiNodes; ++n) {
        *(nodeMappingPtr + n) = n;
    }

    auto* stateIdxPtr = reinterpret_cast<uint32_t*>(instanceData.data() + osdiDesc->state_idx_off);
    for (uint32_t s = 0; s < osdiDesc->num_states; ++s) {
        *(stateIdxPtr + s) = s;
    }

    stateHistN.assign(osdiDesc->num_states, 0.0);
    stateHistNc.assign(osdiDesc->num_states, 0.0);
    stateCheckpoint.assign(osdiDesc->num_states, 0.0);

    osdiReady = true;
}

int BsimCmgFet::getExtraRowRequirement() const {
    if (!osdiReady || !osdiDesc) return 0;
    if (osdiDesc->num_nodes > osdiDesc->num_terminals) {
        return static_cast<int>(osdiDesc->num_nodes - osdiDesc->num_terminals);
    }
    return 0;
}

void BsimCmgFet::allocateRows(int& globalExtraRowCounter) {
    int req = getExtraRowRequirement();
    if (req > 0) {
        assignedExtraRowIdx = MnaIndex::fromExtraRow(globalExtraRowCounter);
        globalExtraRowCounter += req;
    }
}

MnaIndex BsimCmgFet::resolveMnaRowIndex(uint32_t osdiNodeIdx, MnaIndex extraIdx) const {
    if (!osdiDesc) return MnaIndex{};
    uint32_t numTerms = osdiDesc->num_terminals;

    // External terminals: 0=Drain, 1=Gate, 2=Source, 3=Body (defaults to Source if 3-pin netlist)
    if (osdiNodeIdx < numTerms) {
        NetNode netNode{0};
        if (osdiNodeIdx < nodes.size()) {
            netNode = nodes.at(osdiNodeIdx);
        } else if (osdiNodeIdx == 3 && nodes.size() >= 3) {
            netNode = nodes.at(2); // Tie Bulk to Source when 3-pin netlist card is used
        }
        return MnaIndex::fromNetNode(netNode);
    }

    // Internal BSIM-CMG nodes (Di, Si, Gi, thermal T, etc.)
    MnaIndex baseExtra = extraIdx.isValid() ? extraIdx : assignedExtraRowIdx;
    if (!baseExtra.isValid()) return MnaIndex{};
    return MnaIndex::fromExtraRow(baseExtra.idx + static_cast<int>(osdiNodeIdx - numTerms));
}

// ============================================================================
// Nonlinear DC / Transient Stamping (OSDI v0.4 + Built-in Phase-1 Fallback)
// ============================================================================

void BsimCmgFet::stamp(MnaMatrixContext& ctx) {
    if (nodes.size() < 3) return;

    if (!osdiReady || !osdiDesc) {
        stampPhase1(ctx);
        return;
    }

    const uint32_t totalNodes = osdiDesc->num_nodes;
    const uint32_t numTerms   = osdiDesc->num_terminals;
    const auto* collapsed     = reinterpret_cast<const bool*>(instanceData.data() + osdiDesc->collapsed_offset);

    std::vector<double> prevSolve(totalNodes, 0.0);
    if (ctx.currentIterV) {
        for (uint32_t n = 0; n < totalNodes; ++n) {
            MnaIndex mnaIdx = resolveMnaRowIndex(n, ctx.assignedExtraIdx);
            prevSolve.at(n) = getSafeExtraValue(ctx.currentIterV, mnaIdx);
        }
    }

    // Prevent singular extra rows for collapsed internal nodes
    for (uint32_t n = numTerms; n < totalNodes; ++n) {
        if (*(collapsed + n)) {
            MnaIndex row = resolveMnaRowIndex(n, ctx.assignedExtraIdx);
            if (row.isValid()) {
                ctx.G.addValue(row.idx, row.idx, 1.0);
            }
        }
    }

    std::vector<double> currentState = (ctx.orderStep == 2) ? stateHistNc : stateHistN;
    const std::vector<double>& prevState = (ctx.orderStep == 2) ? stateHistNc : stateHistN;

    double alpha = 0.0;
    uint32_t simFlags = CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN | CALC_RESIST_LIM_RHS;
    if (ctx.orderStep == 0) {
        simFlags |= ANALYSIS_DC;
    } else {
        simFlags |= ANALYSIS_TRAN | CALC_REACT_RESIDUAL | CALC_REACT_JACOBIAN | CALC_REACT_LIM_RHS;
        if (ctx.orderStep == 1) {
            alpha = 2.0 / std::max(TR_GAMMA * ctx.dt, 1e-18);
        } else {
            double dt2 = std::max((1.0 - TR_GAMMA) * ctx.dt, 1e-18);
            alpha = (2.0 - TR_GAMMA) / ((1.0 - TR_GAMMA) * dt2);
        }
    }

    DefaultOsdiSimParas simEnv;
    OsdiSimInfo simInfo{};
    simInfo.paras      = simEnv.paras;
    simInfo.abstime    = ctx.time;
    simInfo.prev_solve = prevSolve.data();
    simInfo.prev_state = const_cast<double*>(prevState.data());
    simInfo.next_state = currentState.data();
    simInfo.flags      = simFlags;

    uint32_t retFlags = osdiDesc->eval(nullptr, instanceData.data(), modelData.data(), &simInfo);
    if (retFlags & EVAL_RET_FLAG_FATAL) {
        stampPhase1(ctx);
        return;
    }

    // Bind Jacobian slots inside instanceData and load resistive + reactive entries
    std::vector<double> resistJac;
    std::vector<double> reactJac;
    bindOsdiJacobianSlots(osdiDesc, instanceData.data(), resistJac, reactJac);

    osdiDesc->load_jacobian_resist(instanceData.data(), modelData.data());
    if (ctx.orderStep > 0) {
        osdiDesc->load_jacobian_react(instanceData.data(), modelData.data(), alpha);
    }

    for (uint32_t e = 0; e < osdiDesc->num_jacobian_entries; ++e) {
        const OsdiJacobianEntry& entry = *(osdiDesc->jacobian_entries + e);
        MnaIndex row = resolveMnaRowIndex(entry.nodes.node_1, ctx.assignedExtraIdx);
        MnaIndex col = resolveMnaRowIndex(entry.nodes.node_2, ctx.assignedExtraIdx);
        if (row.isValid() && col.isValid()) {
            double gVal = resistJac.at(e) + ((ctx.orderStep > 0) ? reactJac.at(e) : 0.0);
            if (gVal != 0.0) {
                ctx.G.addValue(row.idx, col.idx, gVal);
            }
        }
    }

    // Load KCL Residuals & Linearized Newton-Raphson RHS Vector into ctx.B
    std::vector<double> resistRhs(totalNodes, 0.0);
    std::vector<double> reactRhs(totalNodes, 0.0);
    osdiDesc->load_residual_resist(instanceData.data(), modelData.data(), resistRhs.data());
    if (ctx.orderStep > 0) {
        osdiDesc->load_residual_react(instanceData.data(), modelData.data(), reactRhs.data());
    }

    for (uint32_t n = 0; n < totalNodes; ++n) {
        MnaIndex row = resolveMnaRowIndex(n, ctx.assignedExtraIdx);
        if (row.isValid()) {
            double fNode = resistRhs.at(n) + ((ctx.orderStep > 0) ? reactRhs.at(n) : 0.0);
            addSafeRhs(ctx.B, row, -fNode);
        }
    }

    // Add J * V_k linearization term to RHS
    for (uint32_t e = 0; e < osdiDesc->num_jacobian_entries; ++e) {
        const OsdiJacobianEntry& entry = *(osdiDesc->jacobian_entries + e);
        MnaIndex row = resolveMnaRowIndex(entry.nodes.node_1, ctx.assignedExtraIdx);
        MnaIndex col = resolveMnaRowIndex(entry.nodes.node_2, ctx.assignedExtraIdx);
        if (row.isValid() && col.isValid() && entry.nodes.node_2 < totalNodes) {
            double gVal = resistJac.at(e) + ((ctx.orderStep > 0) ? reactJac.at(e) : 0.0);
            addSafeRhs(ctx.B, row, gVal * prevSolve.at(entry.nodes.node_2));
        }
    }
}

void BsimCmgFet::updateHistory(const MnaMatrixContext& ctx) {
    if (!osdiReady || !osdiDesc || !ctx.currentIterV) return;

    const uint32_t totalNodes = osdiDesc->num_nodes;
    std::vector<double> convergedV(totalNodes, 0.0);
    for (uint32_t n = 0; n < totalNodes; ++n) {
        MnaIndex mnaIdx = resolveMnaRowIndex(n, ctx.assignedExtraIdx);
        convergedV.at(n) = getSafeExtraValue(ctx.currentIterV, mnaIdx);
    }

    // Use separate temporary buffer to avoid aliasing prev_state and next_state
    std::vector<double> nextState(osdiDesc->num_states, 0.0);
    DefaultOsdiSimParas simEnv;
    OsdiSimInfo simInfo{};
    simInfo.paras      = simEnv.paras;
    simInfo.abstime    = ctx.time;
    simInfo.prev_solve = convergedV.data();
    simInfo.prev_state = stateHistN.data();
    simInfo.next_state = nextState.data();
    simInfo.flags      = ((ctx.orderStep == 0) ? ANALYSIS_DC : ANALYSIS_TRAN) | CALC_REACT_RESIDUAL;

    osdiDesc->eval(nullptr, instanceData.data(), modelData.data(), &simInfo);

    if (ctx.orderStep == 1) {
        stateHistNc = std::move(nextState);
    } else {
        stateHistN = std::move(nextState);
    }
}

// ============================================================================
// Small-Signal Multi-Frequency AC Stamping (OSDI v0.4 + Fallback)
// ============================================================================

void BsimCmgFet::stampAC(CsrMatrixComplex& A_c, std::vector<std::complex<double>>& B_c,
                         double omega, const std::vector<double>& dcState) {
    (void)B_c;
    if (nodes.size() < 3) return;

    if (!osdiReady || !osdiDesc) {
        stampACPhase1(A_c, omega, dcState);
        return;
    }

    const uint32_t totalNodes = osdiDesc->num_nodes;
    const uint32_t numTerms   = osdiDesc->num_terminals;
    const auto* collapsed     = reinterpret_cast<const bool*>(instanceData.data() + osdiDesc->collapsed_offset);

    std::vector<double> dcSolve(totalNodes, 0.0);
    for (uint32_t n = 0; n < totalNodes; ++n) {
        MnaIndex mnaIdx = resolveMnaRowIndex(n, assignedExtraRowIdx);
        dcSolve.at(n) = getSafeExtraValue(dcState, mnaIdx);
    }

    for (uint32_t n = numTerms; n < totalNodes; ++n) {
        if (*(collapsed + n)) {
            MnaIndex row = resolveMnaRowIndex(n, assignedExtraRowIdx);
            if (row.isValid()) {
                A_c.addValue(row.idx, row.idx, std::complex<double>(1.0, 0.0));
            }
        }
    }

        std::vector<double> dummyState(osdiDesc->num_states, 0.0);
    DefaultOsdiSimParas simEnv;
    OsdiSimInfo simInfo{};
    simInfo.paras      = simEnv.paras;
    simInfo.abstime    = 0.0;
    simInfo.prev_solve = dcSolve.data();
    simInfo.prev_state = dummyState.data();
    simInfo.next_state = dummyState.data();
    simInfo.flags      = ANALYSIS_AC | CALC_RESIST_JACOBIAN | CALC_REACT_JACOBIAN;

    uint32_t retFlags = osdiDesc->eval(nullptr, instanceData.data(), modelData.data(), &simInfo);
    if (retFlags & EVAL_RET_FLAG_FATAL) {
        stampACPhase1(A_c, omega, dcState);
        return;
    }

    std::vector<double> resistJac;
    std::vector<double> reactJac;
    bindOsdiJacobianSlots(osdiDesc, instanceData.data(), resistJac, reactJac);

    osdiDesc->load_jacobian_resist(instanceData.data(), modelData.data());
    osdiDesc->load_jacobian_react(instanceData.data(), modelData.data(), 1.0);

    for (uint32_t e = 0; e < osdiDesc->num_jacobian_entries; ++e) {
        const OsdiJacobianEntry& entry = *(osdiDesc->jacobian_entries + e);
        MnaIndex row = resolveMnaRowIndex(entry.nodes.node_1, assignedExtraRowIdx);
        MnaIndex col = resolveMnaRowIndex(entry.nodes.node_2, assignedExtraRowIdx);
        if (row.isValid() && col.isValid()) {
            std::complex<double> yVal(resistJac.at(e), omega * reactJac.at(e));
            if (yVal.real() != 0.0 || yVal.imag() != 0.0) {
                A_c.addValue(row.idx, col.idx, yVal);
            }
        }
    }
}

// ============================================================================
// Sparsity Pattern, Telemetry & Checkpointing
// ============================================================================

void BsimCmgFet::getSparsityPattern(std::vector<std::pair<MnaIndex, MnaIndex>>& patterns,
                                    MnaIndex extraIdx) const {
    if (nodes.size() < 3) return;

    if (osdiReady && osdiDesc) {
        const uint32_t totalNodes = osdiDesc->num_nodes;
        const uint32_t numTerms   = osdiDesc->num_terminals;
        const auto* collapsed     = reinterpret_cast<const bool*>(instanceData.data() + osdiDesc->collapsed_offset);

        for (uint32_t n = numTerms; n < totalNodes; ++n) {
            if (*(collapsed + n)) {
                MnaIndex row = resolveMnaRowIndex(n, extraIdx);
                if (row.isValid()) {
                    patterns.push_back({row, row});
                }
            }
        }

        for (uint32_t e = 0; e < osdiDesc->num_jacobian_entries; ++e) {
            const OsdiJacobianEntry& entry = *(osdiDesc->jacobian_entries + e);
            MnaIndex row = resolveMnaRowIndex(entry.nodes.node_1, extraIdx);
            MnaIndex col = resolveMnaRowIndex(entry.nodes.node_2, extraIdx);
            if (row.isValid() && col.isValid()) {
                patterns.push_back({row, col});
            }
        }
        return;
    }

    // Fallback 3x3 terminal sparsity stencil using type-safe MnaIndex
    for (size_t i : {0, 1, 2}) {
        MnaIndex mi = MnaIndex::fromNetNode(nodes.at(i));
        if (!mi.isValid()) continue;

        for (size_t j : {0, 1, 2}) {
            MnaIndex mj = MnaIndex::fromNetNode(nodes.at(j));
            if (mj.isValid()) {
                patterns.push_back({mi, mj});
            }
        }
    }
}

double BsimCmgFet::calculateVoltageDrop(const std::vector<double>& state) const {
    if (nodes.size() < 3) return 0.0;
    return getSafeVoltage(state, nodes.at(0)) - getSafeVoltage(state, nodes.at(2));
}

double BsimCmgFet::calculateCurrent(const SimulationContext& simCtx, MnaIndex extraIdx) const {
    if (simCtx.state.empty() || nodes.size() < 3) return 0.0;

    if (osdiReady && osdiDesc) {
        const uint32_t totalNodes = osdiDesc->num_nodes;
        std::vector<double> solveV(totalNodes, 0.0);
        for (uint32_t n = 0; n < totalNodes; ++n) {
            MnaIndex mnaIdx = resolveMnaRowIndex(n, extraIdx);
            solveV.at(n) = getSafeExtraValue(simCtx.state, mnaIdx);
        }

        // Local mutable copies ensure const-correctness and thread safety
        std::vector<uint8_t> instCopy = instanceData;
        std::vector<uint8_t> modCopy  = modelData;
        std::vector<double> dummyState(osdiDesc->num_states, 0.0);

        DefaultOsdiSimParas simEnv;
        OsdiSimInfo simInfo{};
        simInfo.paras      = simEnv.paras;
        simInfo.abstime    = 0.0;
        simInfo.prev_solve = solveV.data();
        simInfo.prev_state = dummyState.data();
        simInfo.next_state = dummyState.data();
        simInfo.flags      = ANALYSIS_DC | CALC_RESIST_RESIDUAL;

        if (!(osdiDesc->eval(nullptr, instCopy.data(), modCopy.data(), &simInfo) & EVAL_RET_FLAG_FATAL)) {
            std::vector<double> resistRhs(totalNodes, 0.0);
            osdiDesc->load_residual_resist(instCopy.data(), modCopy.data(), resistRhs.data());
            return resistRhs.front(); // Terminal 0 is Drain current
        }
    }

    // Built-in Phase-1 fallback current evaluation via NetNode getSafeVoltage helper
    double vD = getSafeVoltage(simCtx.state, nodes.at(0));
    double vG = getSafeVoltage(simCtx.state, nodes.at(1));
    double vS = getSafeVoltage(simCtx.state, nodes.at(2));

    const double pol = params.isNChannel ? 1.0 : -1.0;
    vD *= pol;
    vG *= pol;
    vS *= pol;

    bool reversedCurrent = !params.isNChannel;
    double vds = vD - vS;
    if (vds < 0.0) {
        std::swap(vD, vS);
        vds = -vds;
        reversedCurrent = !reversedCurrent;
    }

    double ids = 0.0, gm = 0.0, gds = GMIN_FLOOR;
    evaluatePhase1Physics(vG - vS, vds, ids, gm, gds);
    return reversedCurrent ? -ids : ids;
}

void BsimCmgFet::resetState() {
    std::fill(stateHistN.begin(), stateHistN.end(), 0.0);
    std::fill(stateHistNc.begin(), stateHistNc.end(), 0.0);
}

void BsimCmgFet::saveCheckpoint() {
    stateCheckpoint = stateHistN;
}

void BsimCmgFet::restoreCheckpoint() {
    stateHistN = stateCheckpoint;
}

void BsimCmgFet::setValue(double newVth0) {
    params.VTH0 = newVth0;
    if (osdiReady) {
        setOsdiParamByName("VTH0", newVth0);
    }
}

// ============================================================================
// Built-In Phase-1 Analytical Fallback Kernel
// ============================================================================

void BsimCmgFet::evaluatePhase1Physics(double vgs, double vds,
                                       double& ids, double& gm, double& gds) const {
    const double Leff = std::max(params.L, 1e-9);
    const double nVt  = std::max(params.NFACTOR, 1.0) * V_THERMAL_300K;

    const double vth = params.VTH0 - params.DVT0 - params.ETA0 * vds;
    const double dvth_dvds = -params.ETA0;

    const double denom_vg = 2.0 * nVt;
    double sp_vg = 0.0, dsp_vg = 0.0;
    safeSoftplus((vgs - vth) / denom_vg, sp_vg, dsp_vg);

    const double vgsteff = denom_vg * sp_vg;
    const double dvgsteff_dvgs = dsp_vg;
    const double dvgsteff_dvds = -dsp_vg * dvth_dvds;

    const double esatL = (2.0 * params.VSAT * Leff) / std::max(params.U0, 1e-6);
    const double sum_vdsat = esatL + vgsteff;
    const double vdsat = (esatL * vgsteff) / sum_vdsat;
    const double dvdsat_dvgsteff = (esatL * esatL) / (sum_vdsat * sum_vdsat);

    const double delta = std::max(params.DELTA, 1e-4);
    const double diff = vdsat - vds - delta;
    const double rad = std::sqrt(diff * diff + 4.0 * delta * vdsat);
    const double vdseff = vdsat - 0.5 * (diff + rad);

    const double dvdseff_dvdsat = 1.0 - 0.5 * (1.0 + (diff + 2.0 * delta) / rad);
    const double dvdseff_dvds   = 0.5 * (1.0 + diff / rad);

    const double dvdseff_total_dvgs = dvdseff_dvdsat * (dvdsat_dvgsteff * dvgsteff_dvgs);
    const double dvdseff_total_dvds = dvdseff_dvdsat * (dvdsat_dvgsteff * dvgsteff_dvds) + dvdseff_dvds;

    const double beta_eff = params.U0 * cox * (weff / Leff);
    const double clm_factor = 1.0 + params.PCLM * (vds - vdseff);
    const double vel_denom = 1.0 + (vdseff / esatL);

    const double num_core = (vgsteff - 0.5 * vdseff) * vdseff;
    const double ids_core = beta_eff * (num_core / vel_denom);
    ids = ids_core * clm_factor;

    const double dids_core_dvgsteff = beta_eff * (vdseff / vel_denom);
    const double dnum_dvdseff = vgsteff - vdseff;
    const double dvel_dvdseff = 1.0 / esatL;
    const double dids_core_dvdseff = beta_eff * (
        (dnum_dvdseff * vel_denom - num_core * dvel_dvdseff) / (vel_denom * vel_denom)
    );

    const double dclm_dvgs = -params.PCLM * dvdseff_total_dvgs;
    const double dclm_dvds =  params.PCLM * (1.0 - dvdseff_total_dvds);

    gm = (dids_core_dvgsteff * dvgsteff_dvgs + dids_core_dvdseff * dvdseff_total_dvgs) * clm_factor
       + ids_core * dclm_dvgs;
    gds = (dids_core_dvgsteff * dvgsteff_dvds + dids_core_dvdseff * dvdseff_total_dvds) * clm_factor
        + ids_core * dclm_dvds;

    if (gm < 0.0) gm = 0.0;
    if (gds < GMIN_FLOOR) gds = GMIN_FLOOR;
}

void BsimCmgFet::stampPhase1(MnaMatrixContext& ctx) {
    if (nodes.size() < 3) return;

    NetNode nD = nodes.at(0);
    NetNode nG = nodes.at(1);
    NetNode nS = nodes.at(2);

    double vD = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nD) : 0.0;
    double vG = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nG) : 0.0;
    double vS = ctx.currentIterV ? getSafeVoltage(ctx.currentIterV, nS) : 0.0;

    const double pol = params.isNChannel ? 1.0 : -1.0;
    vD *= pol;
    vG *= pol;
    vS *= pol;

    double vds = vD - vS;
    if (vds < 0.0) {
        std::swap(nD, nS);
        std::swap(vD, vS);
        vds = -vds;
    }
    const double vgs = vG - vS;

    double ids = 0.0, gm = 0.0, gds = GMIN_FLOOR;
    evaluatePhase1Physics(vgs, vds, ids, gm, gds);

    const double ieq = pol * (ids - gm * vgs - gds * vds);

        const MnaIndex mD = MnaIndex::fromNetNode(nD);
    const MnaIndex mG = MnaIndex::fromNetNode(nG);
    const MnaIndex mS = MnaIndex::fromNetNode(nS);

    if (mD.isValid()) {
        ctx.G.addValue(mD.idx, mD.idx,  gds);
        if (mG.isValid()) ctx.G.addValue(mD.idx, mG.idx,  gm);
        if (mS.isValid()) ctx.G.addValue(mD.idx, mS.idx, -(gm + gds));
        addSafeRhs(ctx.B, mD, -ieq);
    }
    if (mS.isValid()) {
        ctx.G.addValue(mS.idx, mS.idx,  (gm + gds));
        if (mD.isValid()) ctx.G.addValue(mS.idx, mD.idx, -gds);
        if (mG.isValid()) ctx.G.addValue(mS.idx, mG.idx, -gm);
        addSafeRhs(ctx.B, mS, ieq);
    }
}

void BsimCmgFet::stampACPhase1(CsrMatrixComplex& A_c, double omega, const std::vector<double>& dcState) {
    if (nodes.size() < 3) return;

    NetNode nD = nodes.at(0);
    NetNode nG = nodes.at(1);
    NetNode nS = nodes.at(2);

    double vD = getSafeVoltage(dcState, nD);
    double vG = getSafeVoltage(dcState, nG);
    double vS = getSafeVoltage(dcState, nS);

    const double pol = params.isNChannel ? 1.0 : -1.0;
    vD *= pol;
    vG *= pol;
    vS *= pol;

    bool swapped = false;
    double vds = vD - vS;
    if (vds < 0.0) {
        std::swap(vD, vS);
        vds = -vds;
        swapped = true;
    }

    double ids = 0.0, gm = 0.0, gds = GMIN_FLOOR;
    evaluatePhase1Physics(vG - vS, vds, ids, gm, gds);

    const double nfinScale = static_cast<double>(std::max(1, params.NFIN));
    std::complex<double> y_gs(0.0, omega * params.CGS0 * nfinScale);
    std::complex<double> y_gd(0.0, omega * params.CGD0 * nfinScale);
    if (swapped) std::swap(y_gs, y_gd);

    const MnaIndex mD = MnaIndex::fromNetNode(swapped ? nS : nD);
    const MnaIndex mG = MnaIndex::fromNetNode(nG);
    const MnaIndex mS = MnaIndex::fromNetNode(swapped ? nD : nS);

    if (mD.isValid()) {
        A_c.addValue(mD.idx, mD.idx, std::complex<double>(gds, 0.0) + y_gd);
        if (mG.isValid()) A_c.addValue(mD.idx, mG.idx, std::complex<double>(gm, 0.0) - y_gd);
        if (mS.isValid()) A_c.addValue(mD.idx, mS.idx, std::complex<double>(-(gm + gds), 0.0));
    }
    if (mG.isValid()) {
        A_c.addValue(mG.idx, mG.idx, y_gs + y_gd);
        if (mD.isValid()) A_c.addValue(mG.idx, mD.idx, -y_gd);
        if (mS.isValid()) A_c.addValue(mG.idx, mS.idx, -y_gs);
    }
    if (mS.isValid()) {
        A_c.addValue(mS.idx, mS.idx, std::complex<double>(gm + gds, 0.0) + y_gs);
        if (mD.isValid()) A_c.addValue(mS.idx, mD.idx, std::complex<double>(-gds, 0.0));
        if (mG.isValid()) A_c.addValue(mS.idx, mG.idx, std::complex<double>(-gm, 0.0) - y_gs);
    }
}



