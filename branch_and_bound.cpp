#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"


static constexpr double K = 10; // max picked assets 
static constexpr double weight_lower_bound = 0.01; //
static constexpr double weight_upper_bound = 1.0; //

struct Individual{
    vector<int> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
};


using Archive = vector<Individual>;
using Point = tuple<double, double>;

using namespace std;

struct Node {
    int level;                        
    int selectedCount;                 

    std::vector<int> fixedInAssets;    

    std::vector<double> LB;  
    double priority;      
    
    bool prunable = false;

    /** TODO: implementar prioridade dos nós */
    bool operator<(const Node& other) const {
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


/*
    TODO: add to UB in the correct position
*/
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
    root.selectedCount = 0;
    root.fixedInAssets.clear();
    root.priority = std::numeric_limits<double>::infinity();

    return root;
}

Node make_with(Node &previous){
    Node with = {};
    with.level = previous.level + 1;
    with.fixedInAssets.push_back(with.level);
    with.selectedCount = previous.selectedCount+1;
    // TODO:  with.priority = ?
    with.priority = 1;

    return with;
}

Node make_without(Node &previous){
    Node without = {};
    without.level = previous.level + 1;
    without.selectedCount = previous.selectedCount;
    // TODO:  with.priority = ?
    without.priority = 1;

    return without;
}



/*
    N is the set of points used to construct a representation of the UB≺
*/
void calculate_nadir_points(vector<Point> nadirPoints, Archive &UB){
    nadirPoints.push_back({UB[0].risk, -INFINITY});
    for(int i = 1; i<UB.size()-1; i++){
        nadirPoints.push_back({UB[i+1].risk, UB[i].expectedReturn});
    }
    nadirPoints.push_back({INFINITY, UB[UB.size()-1].expectedReturn});
}

void calculate_lambdas(vector<double>& lambda, Archive& UB){
    for(int i=0; i<UB.size()-1; i++){
        auto p1 = UB[i];
        auto p2 = UB[i+1];

        auto v11 = p1.risk;
        auto v12 = - p1.expectedReturn;
        auto v21 = p2.risk;
        auto v22 = -p2.expectedReturn;


        lambda[i] = (v22 - v12) / (v11 - v21 + v22 - v12);
    }
}

bool solve(Individual individual, double lambda, Node node){

}

bool test_pruning(vector<tuple<Individual, double>> LB, Archive UB){

}

void bound(Node &node, Archive& UB){
    vector<double> lambdaList = {};
    lambdaList.reserve(UB.size());
    calculate_lambdas(lambdaList, UB);

    vector<tuple<Individual, double>> LB;

    for(auto lambda : lambdaList){
        Individual individual = {};
        bool found = solve(individual, lambda, node);

        if(!found) continue;
        LB.push_back({individual, lambda});

        if(!is_dominated_by_archive(UB, individual)){
            add_to_archive(UB, individual);
        }
    }

    node.prunable = test_pruning(LB, UB);
}

double branch_and_bound(Archive &UB){

    Node root = make_root_node(NUMBER_OF_ASSETS);

    priority_queue<Node> pq;

    /**
     * UB< is the set of points in the object space that is not dominated by any point in UB. 
     * This set can be defined as UB< = {v ∈ R^p | ∀u ∈ UB, u ⊀ v }
     * as this is a continous region, we need a set N (nadirPoints) that limits UB<. So we can define 
     * UB< as UB< = {v ∈ R^p | ∃w ∈ N, v ≺≺ w }
     */
    vector<Point> nadirPoints = {};
    nadirPoints.reserve(UB.size());
    calculate_nadir_points(nadirPoints, UB);

    vector<double> lambdaList = {};
    lambdaList.reserve(UB.size());
    

    pq.push(root);

    while(!pq.empty()){
        Node current = pq.top();
        pq.pop();

        if(current.fixedInAssets.size() > K){ 
            continue;
        }

        // prune this node
        int remaining = NUMBER_OF_ASSETS - current.level +1;
        if(current.fixedInAssets.size() + remaining < K){ 
            continue;
        }

        Node with = make_with(current);
        bound(with, UB);
        if(!with.prunable){
            pq.push(with);
        }
     
        Node without = make_without(current);
        bound(without, UB);
        if(!without.prunable){
            pq.push(without);
        }
 
    }

    return 0.0;
}


int main(){

    Archive UB = convert_nsgaii_population(run_nsgaII()); // best solutions until now
    sort(UB.begin(), UB.end(), [](Individual a, Individual b){
        if(a.risk != b.risk){
            return a.risk < b.risk;
        }
        return a.expectedReturn > b.expectedReturn;
    });

    portfolioData = PortfolioDataLoader::load_from_file("port5.txt"); // TODO: share this variable with the other files

    branch_and_bound(UB);
    return 0;
}