#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include "gurobi_c++.h"
#include "portfolio_data.hpp"
#include <iomanip>

using namespace std;

struct Individual{
    vector<int> picked;
    vector<double> weights;
    double expectedReturn = 0.0; // maximize
    double risk = 0.0;           // minimize (stddev)
};

static inline bool nearly_equal(double a, double b, double eps = 1e-12) {
    return fabs(a - b) <= eps;
}

bool dominates(const Individual& a, const Individual& b){
    const bool betterReturn =
        (a.expectedReturn > b.expectedReturn) &&
        !nearly_equal(a.expectedReturn, b.expectedReturn);

    const bool betterRisk =
        (a.risk < b.risk) &&
        !nearly_equal(a.risk, b.risk);

    const bool noWorseReturn =
        (a.expectedReturn > b.expectedReturn) ||
        nearly_equal(a.expectedReturn, b.expectedReturn);

    const bool noWorseRisk =
        (a.risk < b.risk) ||
        nearly_equal(a.risk, b.risk);

    return (betterReturn && noWorseRisk)
        || (betterRisk && noWorseReturn);
}


static inline void pareto_insert(std::vector<Individual>& pareto,
                                const Individual& cand,
                                double eps = 1e-12)
{
    // 1) If any existing point dominates cand, discard cand
    for (const auto& p : pareto) {
        if (dominates(p, cand)) {
            return;
        }
    }

    // 2) Remove any points dominated by cand
    pareto.erase(
        std::remove_if(pareto.begin(), pareto.end(),
                       [&](const Individual& p) {
                           return dominates(cand, p);
                       }),
        pareto.end()
    );

    // 3) Insert cand
    pareto.push_back(cand);
}


static inline double compute_expected_return_from_w(const vector<double>& w, const PortfolioData& data) {
    double er = 0.0;
    for (int i = 0; i < data.n; ++i) er += w[i] * data.mean[i];
    return er;
}

static inline double compute_variance_from_w(const vector<double>& w, const PortfolioData& data) {
    double var = 0.0;
    const int N = data.n;
    for (int i = 0; i < N; ++i) {
        const double wi = w[i];
        if (nearly_equal(wi, 0.0)) continue;
        for (int j = 0; j < N; ++j) {
            const double wj = w[j];
            if (nearly_equal(wj, 0.0)) continue;
            var += wi * wj * data.cov_at(i, j);
        }
    }
    if (var < 0.0 && var > -1e-10) var = 0.0;
    return var;
}

static inline Individual extract_individual_solution(
    const vector<GRBVar>& y,
    const vector<GRBVar>& w,
    const PortfolioData& data
) {
    const int N = data.n;
    Individual ind;
    ind.picked.assign(N, 0);
    ind.weights.assign(N, 0.0);

    for (int i = 0; i < N; ++i) {
        const double yi = y[i].get(GRB_DoubleAttr_X);
        const double wi = w[i].get(GRB_DoubleAttr_X);

        ind.picked[i]  = (yi > 0.5) ? 1 : 0;
        ind.weights[i] = (wi < 0.0 && wi > -1e-12) ? 0.0 : wi;
    }

    // Compute objectives (return and stddev)
    ind.expectedReturn = compute_expected_return_from_w(ind.weights, data);
    const double var   = compute_variance_from_w(ind.weights, data);
    ind.risk           = std::sqrt(std::max(0.0, var));
    return ind;
}

static inline GRBQuadExpr build_variance_qexpr(const vector<GRBVar>& w, const PortfolioData& data) {
    GRBQuadExpr q = 0.0;
    const int N = data.n;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const double cov = data.cov_at(i, j);
            if (nearly_equal(cov, 0.0)) continue;
            q += cov * w[i] * w[j];
        }
    }
    return q;
}


static inline bool same_solution_weights(const Individual& a, const Individual& b, double tol = 1e-9) {
    if (a.weights.size() != b.weights.size()) return false;
    for (size_t i = 0; i < a.weights.size(); ++i) {
        if (fabs(a.weights[i] - b.weights[i]) > tol) return false;
    }
    return true;
}

vector<Individual> non_dominated_points(
    const PortfolioData& data,
    int K,
    double deltaReturn = 1e-4,
    int maxPoints = 10000,
    int verbose = 0
) {
    vector<Individual> pareto;
    pareto.reserve(256);

    try {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, verbose ? 1 : 0);
        env.start();

        GRBModel model(env);
        const int N = data.n;

        vector<GRBVar> y(N), w(N);
        for (int i = 0; i < N; ++i) {
            y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + std::to_string(i));
            w[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "w_" + std::to_string(i));
        }
        model.update();

        // sum(w)=1
        {
            GRBLinExpr sumW = 0.0;
            for (int i = 0; i < N; ++i) sumW += w[i];
            model.addConstr(sumW == 1.0, "budget");
        }

        // sum(y) <= K
        {
            GRBLinExpr sumY = 0.0;
            for (int i = 0; i < N; ++i) sumY += y[i];
            model.addConstr(sumY <= (double)K, "cardinality");
        }

        // linking: w_i <= y_i
        for (int i = 0; i < N; ++i) {
            model.addConstr(w[i] <= y[i], "link_" + std::to_string(i));
        }

        GRBLinExpr expectedReturnExpr = 0.0;
        for (int i = 0; i < N; ++i) expectedReturnExpr += data.mean[i] * w[i];

        GRBQuadExpr varianceExpr = build_variance_qexpr(w, data);

        GRBConstr fixRetLo; bool hasFixRetLo = false;
        GRBConstr fixRetHi; bool hasFixRetHi = false;

        GRBConstr retUBConstr; bool hasRetUBConstr = false;
        double currentRetUB = std::numeric_limits<double>::infinity();

        const double tolRetFix = 1e-8;
        const double tolDupObj = 1e-9;   // for comparing (return,var)

        for (int iter = 0; iter < maxPoints; ++iter) {
            if (hasFixRetLo) { model.remove(fixRetLo); hasFixRetLo = false; }
            if (hasFixRetHi) { model.remove(fixRetHi); hasFixRetHi = false; }
            model.update();

            // Phase 1: maximize return
            model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
            model.setObjective(expectedReturnExpr);
            model.optimize();

            int status = model.get(GRB_IntAttr_Status);
            if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD) break;
            if (status != GRB_OPTIMAL) break;

            const double bestReturn = expectedReturnExpr.getValue();

            // Phase 2: minimize variance at that bestReturn
            fixRetLo = model.addConstr(expectedReturnExpr >= bestReturn - tolRetFix,
                                       "fixRetLo_" + std::to_string(iter));
            fixRetHi = model.addConstr(expectedReturnExpr <= bestReturn + tolRetFix,
                                       "fixRetHi_" + std::to_string(iter));
            hasFixRetLo = hasFixRetHi = true;
            model.update();

            model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model.setObjective(varianceExpr);
            model.optimize();

            status = model.get(GRB_IntAttr_Status);
            if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD) break;
            if (status != GRB_OPTIMAL) break;

            const double bestVar = varianceExpr.getValue();
            Individual ind = extract_individual_solution(y, w, data);

             const double newUB = ind.expectedReturn - deltaReturn; // use saved return

            if (!pareto.empty()) {
                const auto& prev = pareto.back();
                if (std::fabs(ind.expectedReturn - prev.expectedReturn) <= tolDupObj &&
                    std::fabs(ind.risk*ind.risk - prev.risk*prev.risk) <= tolDupObj) {
                    if (verbose) std::cerr << "Duplicate (obj space). Increase deltaReturn.\n";
                    continue;
                }
            }

            pareto_insert(pareto, ind);

            if (verbose) {
                std::cout << "Point #" << pareto.size()
                          << "  Return=" << ind.expectedReturn
                          << "  Std=" << ind.risk
                          << "  Var=" << bestVar
                          << "  UB(prev)=" << currentRetUB
                          << "\n";
            }

            // Cut return for next iteration:
           
            if (newUB <= 0.0) break;

            if (hasRetUBConstr) {
                model.remove(retUBConstr);
                hasRetUBConstr = false;
                model.update();
            }

            currentRetUB = newUB;
            retUBConstr = model.addConstr(expectedReturnExpr <= currentRetUB,
                                          "cutRet_" + std::to_string(iter));
            hasRetUBConstr = true;
            model.update();
        }

        return pareto;
    } catch (...) {
        throw;
    }
}

int main() {
    PortfolioData data = PortfolioDataLoader::load_from_file("port1.txt");

    int K = 10;

    double deltaReturn = 1e-4;

    auto pareto = non_dominated_points(data, K, deltaReturn, 1000);

    cout << "Found " << pareto.size() << " points\n";
    for (const auto& ind : pareto) {
        cout << "Return=" << ind.expectedReturn << "  Risk(std)=" << ind.risk << "\n";
    }
}