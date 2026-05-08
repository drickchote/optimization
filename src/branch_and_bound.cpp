#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <queue>
#include <string>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"
#include "gurobi_c++.h"
#include <cmath>
#include <map>
#include <chrono>


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
int shit_solutions = 0;


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

double calculate_lambda(const Individual& ind1, const Individual& ind2) {
    double v11 = ind1.risk;
    double v12 = -ind1.expectedReturn;
    double v21 = ind2.risk;
    double v22 = -ind2.expectedReturn;

    double numerator = v22 - v12;
    double denominator = v11 - v21 + v22 - v12;

    if (!std::isfinite(v11) || !std::isfinite(v12) ||
        !std::isfinite(v21) || !std::isfinite(v22)) {
                cout << "!std::isfinite(v11) || !std::isfinite(v12) || !std::isfinite(v21) || !std::isfinite(v22))" << endl;
    }

    if (std::abs(denominator) < 1e-12) {
        cout << "std::abs(denominator) < 1e-12" << endl;
        exit(1);
    }

    double lambda = numerator / denominator;

    if (!std::isfinite(lambda)) {
        cout << "Lambda is infinity" << endl;
        exit(1);

    }

    if (lambda <= 0.0 || lambda >= 1.0) {
        std::cout << "Invalid lambda: " << lambda << std::endl;
        std::cout << "p1: risk=" << ind1.risk
                  << ", return=" << ind1.expectedReturn << std::endl;
        std::cout << "p2: risk=" << ind2.risk
                  << ", return=" << ind2.expectedReturn << std::endl;
        std::cout << "numerator=" << numerator
                  << ", denominator=" << denominator << std::endl;

        exit(1);
    }

    return lambda;
}

void calculate_lambdas(vector<double>& lambda, Archive& UB){
    if (UB.size() < 2) {
        return;
    }
    
    lambda.reserve(UB.size());
    lambda.clear();

    for(int i=0; i<UB.size()-1; i++){
        auto p1 = UB[i];
        auto p2 = UB[i+1];

        double result = calculate_lambda(p1, p2);

        if (result > 0.0 && result < 1.0){
            lambda.push_back(result);
        }
    }
}

/*
    N is the set of points used to construct a representation of the UB≺
*/
void calculate_nadir_points(vector<Point> &nadirPoints, Archive &UB){
    if (UB.size() < 2) {
        return;
    }

    nadirPoints.reserve(UB.size());
    nadirPoints.clear();
    if (UB.size() < 2) return;
    for(int i = 0; i<UB.size()-1; i++){
        nadirPoints.push_back({UB[i+1].risk, -UB[i].expectedReturn});
    }
}

void print_individual(Individual individual){
        cout << setprecision(17) << individual.risk << " " << individual.expectedReturn;

    cout << endl;
}

/**
 * After adding a new individual to UB we need update the lambdas in the following way:
 * Before: p1 λ0 p2 | After: p1 λ1 newP λ2 p2
 * 
 */
void adjust_lambda_after_adding(Archive &UB, vector <double>&lambdaList, int pointPosition, Individual newIndividual){

    /**
     * For calculate lambda we need at last 2 points
     */
    if (UB.size() <= 1) {
        lambdaList.clear();
        return;
    }

    /**
     * It it's first point in UB we add only one lambda between this new point and the next
     */
    if (pointPosition == 0) {
        lambdaList.insert(lambdaList.begin(), calculate_lambda(newIndividual, UB[pointPosition+1]));
        return;
    }


    /**
     * If it's the last point in UB we add only one lambda between this new point and the previous
     */
    if (pointPosition == UB.size()-1) {
        lambdaList.push_back(calculate_lambda(UB[pointPosition -1], newIndividual));
        return;
    }

    auto previousPoint = UB[pointPosition -1];
    auto nextPoint = UB[pointPosition+1];

    /**
     * If the new point is in the middle of UB we need to replace 1 lambda (the one between newIndividual and the previous )
     * and add a new one between the newIndividual and the next point.
     */
    double lambdaBefore = calculate_lambda(previousPoint, newIndividual);
    double lambdaAfter = calculate_lambda(newIndividual, nextPoint);

    lambdaList[pointPosition] = lambdaBefore; // Replace the older lambda for this position
    lambdaList.insert(lambdaList.begin() + pointPosition, lambdaAfter);
}

/**
 * After removing an individual from UB we need update the lambdas in the following way:
 * Before: p1 λ1 p2 λ2 p3 | After: p1 λ0 p3
 * 
 */
void adjust_lambda_after_removal(Archive &UB, vector <double>&lambdaList, int pointPosition){

    if (UB.empty()) {
        lambdaList.clear();
        return;
    }

    if (pointPosition == 0) {
        if (!lambdaList.empty()) {
            lambdaList.erase(lambdaList.begin());
        }
        return;
    }

     // Remotion of the last lambda
    if (pointPosition == UB.size() - 1) {
        lambdaList.pop_back();
        return;
    } 

    lambdaList.erase(lambdaList.begin() + pointPosition);

    auto previousPoint = UB[pointPosition -1];
    auto nextPoint = UB[pointPosition]; // The next point has the same position as the individual that was removed

    double newLambda = calculate_lambda(previousPoint, nextPoint);
    
    lambdaList[pointPosition-1] = newLambda; // Replace the older lambda for this position
}

void adjust_nadir_after_adding(Archive &UB, vector <Point>&nadirList, int pointPosition, Individual newIndividual){

    if (UB.size() <= 1) {
        nadirList.clear();
        return;
    }

    if (pointPosition == 0) {
        nadirList.insert(nadirList.begin(), {UB[pointPosition+1].risk, -UB[pointPosition].expectedReturn});
        return;
    }

    if (pointPosition == UB.size()-1) {
        nadirList.push_back({newIndividual.risk, -UB[pointPosition].expectedReturn});
        return;
    }


    auto previousPoint = UB[pointPosition -1];
    auto nextPoint = UB[pointPosition+1];

    Point nadirBefore = {newIndividual.risk, -previousPoint.expectedReturn};
    Point nadirAfter = {nextPoint.risk, -newIndividual.expectedReturn};


    nadirList[pointPosition] = nadirBefore; // Replace the older lambda for this position
    nadirList.insert(nadirList.begin() + pointPosition, nadirAfter);
}

void adjust_nadir_after_removal(Archive &UB, vector <Point>&nadirList, int pointPosition){
    
    if (UB.empty()) {
        nadirList.clear();
        return;
    }

    if (pointPosition == 0) {
        if (!nadirList.empty()) {
            nadirList.erase(nadirList.begin());
        }
        return;
    }

     // Remotion of the last lambda
    if (pointPosition == UB.size() - 1) {
        nadirList.pop_back();
        return;
    } 

    nadirList.erase(nadirList.begin() + pointPosition);

    auto previousPoint = UB[pointPosition -1];
    auto nextPoint = UB[pointPosition]; // The next point has the same position as the individual that was removed

    Point newNadir = {nextPoint.risk, -previousPoint.expectedReturn};
    
    nadirList[pointPosition-1] = newNadir; // Replace the older nadir for this position
}



void add_to_archive(Archive& UB, const Individual& individual, vector<double> &lambdaList, vector<Point> &nadirPoints) {

    if(UB.size() > 5000){ // Temp TODO - Remove
        return;
    }

    for (auto it = UB.begin(); it != UB.end(); ) {
        if (dominates(individual, *it)) {
            int pointPosition = it - UB.begin();
            it = UB.erase(it);  
            adjust_lambda_after_removal(UB, lambdaList, pointPosition);
            adjust_nadir_after_removal(UB, nadirPoints, pointPosition);

        } else {
            ++it;
        }
    }

    auto it = UB.begin();
    while(it != UB.end() && it->risk < individual.risk){
        ++it;
    }

    int pointPosition = it - UB.begin();
    UB.insert(it, individual);
    adjust_lambda_after_adding(UB, lambdaList, pointPosition, individual);
    adjust_nadir_after_adding(UB, nadirPoints, pointPosition, individual);
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

    return variance;
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

bool dominated_or_existing_solution(Archive& archive, const Individual& ind){
    for (const auto& ub : archive){
        if (dominates(ub, ind) || (nearly_equal(ub.expectedReturn, ind.expectedReturn))){
            if(nearly_equal(ub.risk, ind.risk)){
                // cout << "same shit solution: ";
                shit_solutions++;
                // cout << setprecision(16) << ind.risk<< " " << ind.expectedReturn << " times:" << shit_solutions << endl;;
            }
            return true;

        }
    }
    return false;
}

Node make_root_node(int nAssets) {
    Node root;

    root.level = 0;
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
    with.fixedInAssets.push_back(1);

    with.selectedCount = previous.selectedCount+1;
    with.priority = with.level;

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
    without.fixedInAssets.push_back(0);
    without.priority = without.level;

    return without;
}


bool test_pruning(vector<WeightedBound>& LB, vector<Point>& N) {
    if(LB.empty()){
        return false;
    }

    for (auto& w : N) {
        double risk = get<0>(w);
        double negReturn  = get<1>(w);

        double hLambda = INFINITY;

        for (auto& lb : LB) {
            double lambda = lb.lambda;

            double scalar = lambda * risk
                          + (1 - lambda) * negReturn;

            double value = scalar - lb.lbValue;

            hLambda = min(hLambda, value);
        }

        if (hLambda > 0) {
            return false;
        }
    }

    cout << "Pruning Node" << endl;

    return true;
}

class PortfolioSolver{
    private: 
        GRBEnv env;
        unique_ptr<GRBModel> model;
        vector<GRBVar> w;
        vector<GRBVar> y;
        vector<GRBConstr> fixedConstrs;
public:
    PortfolioSolver() : env(true) {
        try{
            
            env.set(GRB_IntParam_OutputFlag, 0);
            env.start();
            model = std::make_unique<GRBModel>(env);
    
            w.resize(NUMBER_OF_ASSETS);
            y.resize(NUMBER_OF_ASSETS);
    
            for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
                w[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                y[i] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
    
            for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
                model->addConstr(w[i] <= y[i]);
                model->addConstr(w[i] >= weight_lower_bound * y[i]);
            }
    
            GRBLinExpr sumWeights = 0.0;
            GRBLinExpr pickedNumber = 0.0;
    
            for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
                sumWeights += w[i];
                pickedNumber += y[i];
            }
    
            model->addConstr(sumWeights == 1.0);
            model->addConstr(pickedNumber <= K);
    
            model->update();
        } catch (GRBException e){
            cout << e.getMessage() << endl;
            exit(1);
        }
    }

    bool solveNodeLambda(Individual& individual, const Node& node, double lambda) {
        try{
            clearFixedConstraints();
            applyFixedConstraints(node);

            setObjective(lambda);
    
            model->optimize();
            if (model->get(GRB_IntAttr_SolCount) == 0) {
                return false;
            }

            individual.weights.clear();
            individual.picked.clear();
    
            for (int i = 0; i < NUMBER_OF_ASSETS; i++) {
                double wi = w[i].get(GRB_DoubleAttr_X);
                double yi = y[i].get(GRB_DoubleAttr_X);

                individual.weights.push_back(wi);
                individual.picked.push_back(yi > 0.5 ? 1 : 0);
            }
        } catch(GRBException e){
            cout << "Gurobi:" << e.getMessage() << endl;
            exit(1);
        }


        individual.expectedReturn = calculate_expected_return(individual);
        individual.risk = calculate_risk(individual);

        return true;
    }

private:
    void clearFixedConstraints() {
        for (auto& c : fixedConstrs) {
            model->remove(c);
        }
        fixedConstrs.clear();
        model->update();
    }

    void applyFixedConstraints(const Node& node) {
        for (int i = 0; i < node.fixedInAssets.size(); i++) {
            if (node.fixedInAssets[i] == 0) {
                fixedConstrs.push_back(model->addConstr(y[i] == 0));
            } else {
                fixedConstrs.push_back(model->addConstr(y[i] == 1));
            }
        }
        model->update();
    }

    void setObjective(double lambda) {
        GRBQuadExpr risk = calculate_risk(w.data());
        GRBLinExpr expectedReturn = calculate_expected_return(w.data());
        double expr;
        if(lambda == 0){
            model->setObjective(
                -expectedReturn,
                GRB_MINIMIZE
            );
        } else if (lambda == 1){
            model->setObjective(
                risk,
                GRB_MINIMIZE
            );
        } else {
            model->setObjective(
                lambda * risk - (1.0 - lambda) * expectedReturn,
                GRB_MINIMIZE
            );
        }
    }
};


bool bound(Node &node, Archive& UB,  vector<double> &lambdaList, vector<Point> &nadirPoints, PortfolioSolver &solver){
    vector<WeightedBound> LB;
    
    vector<Individual> candidates;
    
    for(auto lambda : lambdaList){
        if (!std::isfinite(lambda)) {
            cout << "Invalid lambda before solve: " << lambda << endl;
            exit(1);
        }
        Individual individual = {};
        bool found = false;
        try {
            found = solver.solveNodeLambda(individual, node, lambda);
        } catch(GRBException e){
            cout << e.getMessage() << endl;
            exit(1);
        }
    
        if(!found) continue;
        double lbValue = lambda * individual.risk - (1 - lambda) * individual.expectedReturn;
        LB.push_back(WeightedBound({lambda, lbValue, individual}));

        candidates.push_back(individual);
    }

    for(auto individual : candidates){
        if(!dominated_or_existing_solution(UB, individual)){
            add_to_archive(UB, individual, lambdaList, nadirPoints);
        } 
    }

    node.prunable = test_pruning(LB, nadirPoints);
}

static constexpr const char* CKPT_MAGIC = "BB_CKPT_V1";

static bool read_individual(std::istream& in, Individual& ind) {
    int pc, wc;
    in >> pc;
    ind.picked.resize(pc);
    for (int k = 0; k < pc; ++k) in >> ind.picked[k];
    in >> wc;
    ind.weights.resize(wc);
    for (int k = 0; k < wc; ++k) in >> ind.weights[k];
    in >> ind.expectedReturn >> ind.risk;
    return static_cast<bool>(in);
}

static void write_individual(std::ostream& out, const Individual& ind) {
    out << ind.picked.size();
    for (int v : ind.picked) out << ' ' << v;
    out << '\n' << ind.weights.size();
    for (double w : ind.weights) out << ' ' << std::setprecision(17) << w;
    out << '\n' << std::setprecision(17) << ind.expectedReturn << ' ' << ind.risk << '\n';
}

static bool read_node(std::istream& in, Node& node) {
    in >> node.level >> node.selectedCount;
    int fc;
    in >> fc;
    node.fixedInAssets.resize(fc);
    for (int k = 0; k < fc; ++k) in >> node.fixedInAssets[k];
    int lbc;
    in >> lbc;
    node.LB.resize(lbc);
    for (int k = 0; k < lbc; ++k) in >> node.LB[k];
    int pr;
    in >> node.priority >> pr;
    node.prunable = (pr != 0);
    return static_cast<bool>(in);
}

static void write_node(std::ostream& out, const Node& node) {
    out << node.level << ' ' << node.selectedCount << '\n';
    out << node.fixedInAssets.size();
    for (int a : node.fixedInAssets) out << ' ' << a;
    out << '\n' << node.LB.size();
    for (double v : node.LB) out << ' ' << std::setprecision(17) << v;
    out << '\n' << std::setprecision(17) << node.priority << ' ' << (node.prunable ? 1 : 0) << '\n';
}

static bool load_checkpoint(const std::string& path, Archive& UB, std::priority_queue<Node>& pq, int& i) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Cannot open checkpoint for reading: " << path << std::endl;
        return false;
    }
    std::string magic;
    in >> magic;
    if (magic != CKPT_MAGIC) {
        std::cerr << "Invalid checkpoint header.\n";
        return false;
    }
    int nAssets = 0;
    in >> nAssets;
    if (nAssets != portfolioData.n) {
        std::cerr << "Checkpoint n_assets (" << nAssets << ") != portfolioData.n (" << portfolioData.n << ").\n";
        return false;
    }
    int saved_i = 0;
    in >> saved_i;
    i = saved_i - 1;
    int ubCount = 0;
    in >> ubCount;
    UB.clear();
    UB.reserve(ubCount);
    for (int u = 0; u < ubCount; ++u) {
        Individual ind{};
        if (!read_individual(in, ind)) return false;
        UB.push_back(std::move(ind));
    }
    int pqCount = 0;
    in >> pqCount;
    while (!pq.empty()) pq.pop();
    for (int p = 0; p < pqCount; ++p) {
        Node node{};
        if (!read_node(in, node)) return false;
        pq.push(std::move(node));
    }
    return static_cast<bool>(in);
}

static void save_checkpoint(const std::string& path, const Archive& UB, const std::priority_queue<Node>& pq, int i) {
    std::priority_queue<Node> copy = pq;
    std::vector<Node> nodes;
    nodes.reserve(copy.size());
    while (!copy.empty()) {
        nodes.push_back(copy.top());
        copy.pop();
    }
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Cannot open checkpoint for writing: " << path << std::endl;
        return;
    }
    out << std::setprecision(17);
    out << CKPT_MAGIC << '\n';
    out << portfolioData.n << '\n';
    out << i << '\n';
    out << static_cast<int>(UB.size()) << '\n';
    for (const Individual& ind : UB) write_individual(out, ind);
    out << static_cast<int>(nodes.size()) << '\n';
    for (const Node& node : nodes) write_node(out, node);
}

double branch_and_bound(Archive &UB, const std::string& checkpoint_in, const std::string& checkpoint_out){
    priority_queue<Node> pq;
    int i = 0;

    if (!checkpoint_in.empty()) {
        if (!load_checkpoint(checkpoint_in, UB, pq, i)) {
            std::cerr << "Aborting: checkpoint load failed.\n";
            return -1.0;
        }
    } else {
        Node root = make_root_node(NUMBER_OF_ASSETS);
        pq.push(root);
    }

    /**
     * UB< is the set of points in the object space that is not dominated by any point in UB. 
     * This set can be defined as UB< = {v ∈ R^p | ∀u ∈ UB, u ⊀ v }
     * as this is a continous region, we need a set N (nadirPoints) that limits UB<. So we can define 
     * UB< as UB< = {v ∈ R^p | ∃w ∈ N, v ≺≺ w }
     */
    vector<Point> nadirPoints = {};
    vector<double> lambdaList = {};

    unsigned long long numberOfNodes = (1ULL << (NUMBER_OF_ASSETS+1)) - 1;
    int treeHeight = NUMBER_OF_ASSETS + 1;
    unsigned long long missingNodes = numberOfNodes;

    calculate_nadir_points(nadirPoints, UB);
    calculate_lambdas(lambdaList, UB);

    PortfolioSolver solver;

    const auto time_start = std::chrono::steady_clock::now();
    while(!pq.empty()){
        cout << "UB size " << UB.size() << endl; 
        cout << "Lambda size" << lambdaList.size() << endl; 
        cout << "nadir size" << nadirPoints.size() << endl; 

        missingNodes--;
        cout << "Missing Nodes: "<< missingNodes << endl;
        i++; 
        if(i % 10 == 0){
            const auto elapsed = std::chrono::steady_clock::now() - time_start;
            const double minutes =
                std::chrono::duration<double>(elapsed).count() / 60.0;
            cout << "minutes since start: " << fixed << setprecision(2) << minutes << endl;
            cout << "i value: " << i << endl;
            cout << "queue size " << pq.size() << endl;
            if (!checkpoint_out.empty()) {
                save_checkpoint(checkpoint_out, UB, pq, i);
            }
        }
        Node current = pq.top();
        pq.pop();

        
        if(current.selectedCount >= K || current.level == NUMBER_OF_ASSETS){ 
            continue;
        }

        Node with = make_with(current);
        bound(with, UB, lambdaList,  nadirPoints, solver);
        if(!with.prunable){
            pq.push(with);
        } else {
            int remotion = (1 << (treeHeight-with.level)) - 1;
            missingNodes -= remotion;
            cout << "A node in the level "<< with.level << " was pruned. It removed " << remotion << " nodes" << endl;
            cout << "Missing Nodes: " << missingNodes << endl;
        }
     
        Node without = make_without(current);
        bound(without, UB, lambdaList, nadirPoints, solver);
        if(!without.prunable){
            pq.push(without);
        } else {
            int remotion = (1 << (treeHeight-with.level)) - 1;
            missingNodes -= remotion;
            cout << "A node in the level "<< with.level << " was pruned. It removed " << remotion << " nodes" << endl;
            cout << "Missing Nodes: " << missingNodes << endl;
        }
 
    }
    cout << "finished BB" << endl;
    return 0.0;
}




int main(int argc, char** argv) {
    if (argc != 1 && argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [checkpoint_in checkpoint_out]\n";
        std::cerr << "  checkpoint_in: file to resume from, or '-' to start from NSGA-II + B&B.\n";
        std::cerr << "  checkpoint_out: file written every 10000 iterations, or '-' to disable.\n";
        return 1;
    }

    std::string checkpoint_in;
    std::string checkpoint_out;
    if (argc == 3) {
        checkpoint_in = argv[1];
        checkpoint_out = argv[2];
        if (checkpoint_in == "-") checkpoint_in.clear();
        if (checkpoint_out == "-") checkpoint_out.clear();
    }

    Archive UB = {};

    if (!checkpoint_in.empty()) {
        portfolioData = PortfolioDataLoader::load_from_file("../inputs/port1.txt");
        NUMBER_OF_ASSETS = portfolioData.n;
    } else {
        convert_nsgaii_population(UB, run_nsgaII());
        sort(UB.begin(), UB.end(), [](Individual a, Individual b) {
            if (a.risk != b.risk) {
                return a.risk < b.risk;
            }
            return a.expectedReturn > b.expectedReturn;
        });
        portfolioData = PortfolioDataLoader::load_from_file("../inputs/port1.txt");
        NUMBER_OF_ASSETS = portfolioData.n;
    }

    if (branch_and_bound(UB, checkpoint_in, checkpoint_out) < 0.0) {
        return 1;
    }
    return 0;
}