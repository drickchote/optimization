#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"


static constexpr double K = 10; // max picked assets 
static constexpr double LB = 0.01; //
static constexpr double UB = 1.0; // 

struct Individual{
    vector<int> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
};


using Archive = vector<Individual>;

using namespace std;

struct Node {
    int level;                        
    std::vector<int8_t> fix;            // -1 = free, 0 = excluded, 1 = included
    int selectedCount;                 

    std::vector<int> fixedInAssets;     // índices com fix[i] == 1

    std::vector<double> LB_lambda;      // LB_λ(N) para cada λ ∈ Λ
    double priority;                    // escalar heurístico p/ fila

    bool operator>(const Node& other) const {
        return priority > other.priority;
    }
};

double calculate_expected_return(const Individual& ind) {
    const int N = portfolioData.n;
    double er = 0.0;

    for (int i = 0; i < N; ++i) {
        er += ind.picked[i] * ind.weights[i] * portfolioData.mean[i];
    }

    return er;
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

void add_to_archive(Archive& UB, const Individual& individual) {
    for (auto it = UB.begin(); it != UB.end(); ) {
        if (dominates(individual, *it)) {
            it = UB.erase(it);  
        } else {
            ++it;
        }
    }

    UB.push_back(individual);
}


double calculate_risk(const Individual& ind) {
    const int N = portfolioData.n;
    double variance = 0.0;

    for (int i = 0; i < N; ++i) {
        if (ind.picked[i] == 0) continue;

        const double wi = ind.weights[i];
        const double si = portfolioData.stdev[i];

        for (int j = 0; j < N; ++j) {
            if (ind.picked[j] == 0) continue;

            const double wj = ind.weights[j];
            const double sj = portfolioData.stdev[j];

            variance += wi * wj * si * sj * portfolioData.corr_at(i, j);
        }
    }

    if (variance < 0.0 && variance > -1e-12)
        variance = 0.0;

    return sqrt(variance);
}


Archive convert_nsgaii_population(NSGAII_Population population){
    Archive archive = {};
    archive.reserve(population.size());

    for(NSGAII_Individual nsgaII_individual : population){
        Individual individual = {};
        individual.expectedReturn = nsgaII_individual.expectedReturn;
        individual.risk = nsgaII_individual.risk;
        individual.picked = nsgaII_individual.picked;
        individual.weights = nsgaII_individual.weights;
    }

    return archive;
}

bool is_dominated_by_archive(const Archive& UB, Individual individual){
    for (auto ubIndividual : UB){
        if(dominates(ubIndividual, individual)){
            return true;
        }
    }
    return false;
}

Node make_root_node(int nAssets) {
    Node root;

    root.level = -1;

    /*
        Maybe I don't need -1. Only 0 and 1
        If that is the case, this array is not really necessary
    */ 
    root.fix.assign(nAssets, -1);

    root.selectedCount = 0;

    root.fixedInAssets.clear();

    root.LB_lambda.clear();

    root.priority = std::numeric_limits<double>::infinity();

    return root;
}

Node make_with(Node &previous){
    Node with = {};
    with.level = previous.level + 1;
    with.fixedInAssets.push_back(with.level);
    with.selectedCount = previous.selectedCount+1;
    with.fix[with.level] = 1; // is this needed?

    return with;
}

Node make_without(Node &previous){
    Node without = {};
    without.level = previous.level + 1;
    without.selectedCount = previous.selectedCount;
    without.fix[without.level] = 0; // is this needed?

    return without;
}

bool build_feasible_weights_with_bounds(Individual& sol) {
    const int K = static_cast<int>(sol.picked.size());
    if (K <= 0) return false;

    if (K * LB > 1.0 + 1e-12) return false;
    if (K * UB < 1.0 - 1e-12) return false;

    // Start all weights with 0
    sol.weights.assign(sol.picked.size(), 0.0);

    // 1) Same weights
    const double w0 = 1.0 / static_cast<double>(K);

    for (int i = 0; i < K; ++i) {
        sol.weights[i] = w0;
    }

    // 2) apply limits
    double sum = 0.0;
    for (int i = 0; i < K; ++i) {
        if (sol.weights[i] < LB) sol.weights[i] = LB;
        if (sol.weights[i] > UB) sol.weights[i] = UB;
        sum += sol.weights[i];
    }

    // 3) normalize
    double diff = 1.0 - sum;

    const int maxIter = 1000;
    int iter = 0;

    while (std::fabs(diff) > 1e-12 && iter < maxIter) {
        bool adjusted = false;

        for (int i = 0; i < K && std::fabs(diff) > 1e-12; ++i) {
            if (diff > 0.0) {
                // Try to improve
                double slack = UB - sol.weights[i];
                if (slack > 0.0) {
                    double delta = std::min(slack, diff);
                    sol.weights[i] += delta;
                    diff -= delta;
                    adjusted = true;
                }
            } else {
                // Try to reduce
                double slack = sol.weights[i] - LB;
                if (slack > 0.0) {
                    double delta = std::min(slack, -diff);
                    sol.weights[i] -= delta;
                    diff += delta;
                    adjusted = true;
                }
            }
        }

        if (!adjusted) break;
        ++iter;
    }

    double finalSum = 0.0;
    for (int i = 0; i < K; ++i) {
        if (sol.weights[i] < LB - 1e-10) return false;
        if (sol.weights[i] > UB + 1e-10) return false;
        finalSum += sol.weights[i];
    }

    if (std::fabs(finalSum - 1.0) > 1e-8) return false;

    return true;
}


bool try_get_feasible_solution_from_node(
    const Node& N,
    Individual& outSolution
) {
    if(N.selectedCount != K) return false;
    
    const int nAssets = static_cast<int>(N.fix.size());
    if (nAssets <= 0) return false;

    outSolution.picked.reserve(nAssets);
    outSolution.picked = std::move(N.fixedInAssets);
  

    outSolution.weights.reserve(nAssets);


    if (!build_feasible_weights_with_bounds(outSolution)) {
        return false;
    }

    outSolution.expectedReturn = calculate_expected_return(outSolution);
    outSolution.risk  = calculate_risk(outSolution);

    return true;
}


double branch_and_bound(Archive &arquive){

    Node root = make_root_node(NUMBER_OF_ASSETS);

    priority_queue<Node> pq;

    pq.push(root);

    
    while(!pq.empty()){
        Node current = pq.top();
        pq.pop();

        // // prune this node
        // I believe this will never be greater than K
        // if(current.fixedInAssets.size() > K){ 
        //     continue;
        // }

        // prune this node
        int remaining = NUMBER_OF_ASSETS - current.level +1;
        if(current.fixedInAssets.size() + remaining < K){ 
            continue;
        }

        Node with = make_with(current);
        if(with.fixedInAssets.size() == K) {
            Individual newSolution = {};
            bool foundFeasible = try_get_feasible_solution_from_node(with, newSolution);
            if(foundFeasible && !is_dominated_by_archive(arquive, newSolution)){
                add_to_archive(arquive, newSolution);
            }
        } else {
            pq.push(with);
        }


        Node without = make_without(current);
        if(without.fixedInAssets.size() == K){
            Individual newSolution = {};
            bool foundFeasible = try_get_feasible_solution_from_node(without, newSolution);
            if(foundFeasible && !is_dominated_by_archive(arquive, newSolution)){
                add_to_archive(arquive, newSolution);
            }
        } else {
            pq.push(without);
        }
 
    }

}










int main(){
    Archive archive = convert_nsgaii_population(run_nsgaII()); // best solutions until now
    branch_and_bound(archive);
    
}