#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"
#include "gurobi_c++.h"
#include <cmath>
#include <map>


static constexpr double K = 10; // max picked assets 
static constexpr double weight_lower_bound = 0.01; //
static constexpr double weight_upper_bound = 1.0; //

struct Individual{
    vector<int> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
};

struct WeightedBound {
    double lambda;
    double lbValue;
    Individual individual;
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

GRBLinExpr calculate_expected_return(const GRBVar weights[]) {
    GRBLinExpr expr = 0.0;

    for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
        expr += weights[i] * portfolioData.mean[i];
    }

    return expr;
}

GRBQuadExpr calculate_risk(const GRBVar weights[]) {
    const int N = portfolioData.n;
    GRBQuadExpr variance = 0.0;


    for (int i = 0; i < N; ++i) {

        const double si = portfolioData.stdev[i];

        for (int j = 0; j < N; ++j) {
            const double sj = portfolioData.stdev[j];

            variance += weights[i] *  weights[j] * si * sj * portfolioData.corr_at(i, j);
        }
    }

    return variance;
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


void convert_nsgaii_population(Archive &UB, NSGAII_Population population){
    UB.reserve(population.size());

    for(NSGAII_Individual nsgaII_individual : population){
        Individual individual = {};
        individual.expectedReturn = nsgaII_individual.expectedReturn;
        individual.risk = nsgaII_individual.risk;
        individual.picked = nsgaII_individual.picked;
        individual.weights = nsgaII_individual.weights;
        UB.push_back(individual);
    }
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
    with.fixedInAssets =  {};

    for(auto asset : previous.fixedInAssets){
        with.fixedInAssets.push_back(asset);
    }
    with.fixedInAssets.push_back(with.level);

    with.selectedCount = previous.selectedCount+1;
    with.priority = 1;

    return with;
}

Node make_without(Node &previous){
    Node without = {};
    without.level = previous.level + 1;
    without.selectedCount = previous.selectedCount;
    without.fixedInAssets = {};
    for(auto asset : previous.fixedInAssets){
        without.fixedInAssets.push_back(asset);
    }
    without.priority = 1;

    return without;
}



/*
    N is the set of points used to construct a representation of the UB≺
*/
void calculate_nadir_points(vector<Point> &nadirPoints, Archive &UB){
    nadirPoints.clear();
    if (UB.size() < 2) return;
    cout << "UB size " << UB.size() << endl; 
    for(int i = 0; i<UB.size()-1; i++){
        nadirPoints.push_back({UB[i+1].risk, -UB[i].expectedReturn});
    }
}

void calculate_lambdas(vector<double>& lambda, Archive& UB){
    for(int i=0; i<UB.size()-1; i++){
        auto p1 = UB[i];
        auto p2 = UB[i+1];

        auto v11 = p1.risk;
        auto v12 = - p1.expectedReturn;
        auto v21 = p2.risk;
        auto v22 = -p2.expectedReturn;

        lambda.push_back((v22 - v12) / (v11 - v21 + v22 - v12));
    }
}


/**
 * Solve min(λf1(x) - (1-λ)f2(x))
 */
bool solve(Individual &individual, double lambda, Node& node){
    try{
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.start();
    
        GRBModel model = GRBModel(env);

        GRBVar* w = model.addVars(NUMBER_OF_ASSETS);

        for(int i = 0; i < NUMBER_OF_ASSETS; i++){
            w[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }

        for(auto assetIndex : node.fixedInAssets){
            model.addConstr(w[assetIndex] >= 0.01);
        }

        GRBLinExpr sumWeights = 0.0;
        for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
            sumWeights += w[i];
        }
        model.addConstr(sumWeights == 1.0);
                
        auto risk = calculate_risk(w);
        auto expectedReturn = calculate_expected_return(w);
        
        
        model.setObjective(lambda * risk - (1-lambda) * expectedReturn, GRB_MINIMIZE);
        model.optimize();

        bool found = false;

        if (model.get(GRB_IntAttr_SolCount) > 0) {
            found = true;

            for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
                individual.weights.push_back(w[i].get(GRB_DoubleAttr_X));
                individual.picked.push_back(individual.weights[i] > 0.0 ? 1 : 0);
            }
            individual.expectedReturn = calculate_expected_return(individual);
            individual.risk = calculate_risk(individual);
        }

        return found;
     } catch (GRBException &e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        exit(1);
    } catch (...) {
        cout << "Exception during optimization" << endl;
        exit(1);
    }
}

bool test_pruning(vector<WeightedBound>& LB, vector<Point>& N) {
    for (auto& w : N) {
        double risk = get<0>(w);
        double negReturn  = get<1>(w);

        double hLambda = INFINITY;

        for (auto& lb : LB) {
            double lambda = lb.lambda;

            double scalar = lambda * risk
                          - (1 - lambda) * negReturn;

            double value = scalar - lb.lbValue;

            hLambda = min(hLambda, value);
        }

        if (hLambda > 0) {
            cout << "Lambda: " << hLambda << endl;
            return false;
        }
    }
    cout << "Pruning Node" << endl;
    return true;
}



bool bound(Node &node, Archive& UB,  vector<double> &lambdaList,  bool shouldCalculateLambdas, vector<Point> &N){
    if(shouldCalculateLambdas){
        lambdaList.clear();
        lambdaList.reserve(UB.size());
        calculate_lambdas(lambdaList, UB);
    }

    vector<WeightedBound> LB;
    bool upperBoundHasChanged = false;
    for(auto lambda : lambdaList){
        Individual individual = {};
        bool found = solve(individual, lambda, node);

    
        if(!found) continue;
        double lbValue = lambda * individual.risk - (1 - lambda) * individual.expectedReturn;
        LB.push_back(WeightedBound({lambda, lbValue, individual}));

        if(!is_dominated_by_archive(UB, individual)){
            add_to_archive(UB, individual);
            upperBoundHasChanged = true;
        } 
    }

    node.prunable = test_pruning(LB, N);
    return upperBoundHasChanged;
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
    vector<double> lambdaList = {};

    pq.push(root);


    bool upperBoundHasChanged = true;

    int i=0;
    while(!pq.empty()){
        i++;
        if(i % 10000 == 0){
            cout << "i value: " << i << endl;
            cout << "queue size " << pq.size() << endl;
        }
        Node current = pq.top();
        pq.pop();
       

        
        if(current.fixedInAssets.size() > K || current.level == NUMBER_OF_ASSETS -1){ 
            continue;
        }

        if(upperBoundHasChanged){
            calculate_nadir_points(nadirPoints, UB);
        }

        Node with = make_with(current);
        upperBoundHasChanged = bound(with, UB, lambdaList, upperBoundHasChanged, nadirPoints);
        if(!with.prunable){
            pq.push(with);
        }
     
        Node without = make_without(current);
        upperBoundHasChanged = bound(without, UB, lambdaList, upperBoundHasChanged, nadirPoints);
        if(!without.prunable){
            pq.push(without);
        }
 
    }
    cout << "finished BB" << endl;
    return 0.0;
}


int main(){

    
    Archive UB = {};
    convert_nsgaii_population(UB, run_nsgaII()); // best solutions until now

    sort(UB.begin(), UB.end(), [](Individual a, Individual b){
        if(a.risk != b.risk){
            return a.risk < b.risk;
        }
        return a.expectedReturn > b.expectedReturn;
    });

    portfolioData = PortfolioDataLoader::load_from_file("port1.txt"); // TODO: share this variable with the other files

    branch_and_bound(UB);
    return 0;
}