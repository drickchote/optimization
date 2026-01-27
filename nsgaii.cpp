#include <iostream>
#include <list>
#include <tuple>
#include <random>
#include "portfolio_data.hpp"

using namespace std;

static int DEBUG = 0;


static constexpr int POP_SIZE = 40;
static constexpr int GENERATIONS = 1;
static constexpr int NUMBER_OF_ASSETS = 31; // get from the file portx.txt
static constexpr double LB = 0.01; //
static constexpr double UB = 1.0; // 
static constexpr double K = 10; // max picked assets 



mt19937 rng(42);


struct Individual{
    vector<int> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
    int rank;
    int dominatedCount = 0;
    vector<int> dominates;
    double crowdingDistance = 0.0;
};

using Population = vector<Individual>;
using Frontiers = vector<vector<Individual>>; 

static inline int randi(mt19937& rng, int lo, int hi_inclusive) {
    uniform_int_distribution<int> d(lo, hi_inclusive);
    return d(rng);
}

static inline double rand01(mt19937& rng) {
    uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng);
}

void print_individual(Individual individual){
    cout << "______________" << endl;
    for(int i=0; i<NUMBER_OF_ASSETS; i++){
        if(individual.picked[i]){
            cout << i <<" - " << individual.weights[i] << " | ";
        }
    }
    cout << endl;

    cout << "Expected Return: " << individual.expectedReturn << endl;
    cout << "Risk: " << individual.risk << endl;
    cout << "______________" << endl;
}

void print_population(Population &population){
    for(auto i : population){
        print_individual(i);
    }
}

void print_frontiers(Frontiers frontiers, Population population){

    int frontierNumber = 1;
    for(auto frontier : frontiers){
        cout << "Frontier "<< frontierNumber << endl;
        frontierNumber++;
        for(auto individual : frontier){
            print_individual(individual);
        }
        cout << "------------" << endl;
    }
}

double calculate_expected_return(const Individual& ind, const PortfolioData& data) {
    const int N = data.n;
    double er = 0.0;

    for (int i = 0; i < N; ++i) {
        if (ind.picked[i] == 0) continue;
        er += ind.weights[i] * data.mean[i];
    }

    return er;
}

double calculate_risk(const Individual& ind, const PortfolioData& data) {
    const int N = data.n;
    double variance = 0.0;

    for (int i = 0; i < N; ++i) {
        if (ind.picked[i] == 0) continue;

        const double wi = ind.weights[i];
        const double si = data.stdev[i];

        for (int j = 0; j < N; ++j) {
            if (ind.picked[j] == 0) continue;

            const double wj = ind.weights[j];
            const double sj = data.stdev[j];

            variance += wi * wj * si * sj * data.corr_at(i, j);
        }
    }

    if (variance < 0.0 && variance > -1e-12)
        variance = 0.0;

    return sqrt(variance);
}


Individual generate_decision(int size, mt19937& rng) {
    Individual individual;


    uniform_int_distribution<int> binaryDist(0, 1);
    uniform_real_distribution<double> realDist(0.0, 1.0);

    individual.picked.resize(size);
    individual.weights.resize(size, 0.0);

    int pickedCount = 0;
    for (int i = 0; i < size; ++i) {
        individual.picked[i] = binaryDist(rng);
        if (individual.picked[i] == 1)
            pickedCount++;
    }

    // At least one should be selected
    if (pickedCount == 0) {
        int idx = uniform_int_distribution<int>(0, size - 1)(rng);
        individual.picked[idx] = 1;
        pickedCount = 1;
    }

    // Generate weights for the selected assets
    double sumWeights = 0.0;
    for (int i = 0; i < size; ++i) {
        if (individual.picked[i] == 1) {
            individual.weights[i] = realDist(rng);
            sumWeights += individual.weights[i];
        }
    }

    // normalization
    for (int i = 0; i < size; ++i) {
        if (individual.picked[i] == 1) {
            individual.weights[i] /= sumWeights;
        }
    }

    return individual;
}

static inline bool nearly_equal(double a, double b, double eps = 1e-12) {
    return fabs(a - b) <= eps;
}


Population generate_population(int size, PortfolioData data){ // generates a random population
    Population population;
    population.reserve(POP_SIZE);

    for(int i=0; i<size; i++){
        Individual individual = generate_decision(NUMBER_OF_ASSETS, rng);
        individual.expectedReturn = calculate_expected_return(individual, data);
        individual.risk = calculate_risk(individual, data);
        population.push_back(individual);
    }

    return population;
}

Population generate_population(Population population){ // generates a population with genetic operations

}

Population crossover(Individual individual1, Individual individual2) {
    Population children;
    children.reserve(2);

    const double alpha = 0.5; // BLX-α 

    Individual child1, child2;
    child1.picked.resize(NUMBER_OF_ASSETS, 0);
    child2.picked.resize(NUMBER_OF_ASSETS, 0);
    child1.weights.resize(NUMBER_OF_ASSETS, 0.0);
    child2.weights.resize(NUMBER_OF_ASSETS, 0.0);

    // 1) picked: gene-a-gene (uniform crossover in the binary part)
    for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
        if (rand01(rng) < 0.5) {
            child1.picked[i] = individual1.picked[i];
            child2.picked[i] = individual2.picked[i];
        } else {
            child1.picked[i] = individual2.picked[i];
            child2.picked[i] = individual1.picked[i];
        }
    }

    auto blx_gene = [&](double a, double b) -> double {
        const double lo = min(a, b);
        const double hi = max(a, b);
        const double d  = hi - lo;

        const double L = lo - alpha * d;
        const double U = hi + alpha * d;

        uniform_real_distribution<double> dist(L, U);
        return dist(rng);
    };

    auto normalize_child = [&](Individual& c) {
        // zerate weights of non chosed
        double sum = 0.0;
        int countPicked = 0;

        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (c.picked[i] == 0) {
                c.weights[i] = 0.0;
                continue;
            }
            countPicked++;

            // clamp
            if (c.weights[i] < LB) c.weights[i] = LB;
            if (c.weights[i] > UB) c.weights[i] = UB;

            sum += c.weights[i];
        }

        // if none was chosed forces 1
        if (countPicked == 0) {
            int idx = randi(rng, 0, NUMBER_OF_ASSETS - 1);
            c.picked[idx] = 1;
            c.weights[idx] = 1.0;
            // resto já está 0
            return;
        }

        // if the sums still 0 forces 1
        if (nearly_equal(sum, 0.0)) {
            int idx = -1;
            for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
                if (c.picked[i]) { idx = i; break; }
            }
            for (int i = 0; i < NUMBER_OF_ASSETS; ++i) c.weights[i] = 0.0;
            c.weights[idx] = 1.0;
            return;
        }

        // normalize to sum 1
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (c.picked[i] == 1) c.weights[i] /= sum;
        }
    };

    // 2) weights:
    for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
        const double p1 = individual1.weights[i];
        const double p2 = individual2.weights[i];

        child1.weights[i] = blx_gene(p1, p2);
        child2.weights[i] = blx_gene(p1, p2);
    }

    normalize_child(child1);
    normalize_child(child2);

    children.push_back(move(child1));
    children.push_back(move(child2));
    return children;
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
        nearly_equal(a.risk, b.expectedReturn);

    return (betterReturn && noWorseRisk)
        || (betterRisk && noWorseReturn);
}


void evaluate(Population& population){
    const int n = (int)population.size();

    for (int i = 0; i < n; i++){
        population[i].dominatedCount = 0;
        population[i].dominates.clear();
    }

    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            if (dominates(population[i], population[j])) {
                population[i].dominates.push_back(j);
                population[j].dominatedCount++;
            } else if (dominates(population[j], population[i])) {
                population[j].dominates.push_back(i);
                population[i].dominatedCount++;
            }
        }
    }
}

Frontiers non_dominant_sort(Population& population){
    evaluate(population);

    vector<vector<int>> fronts_idx;
    vector<int> current;

    const int n = (int)population.size();
    current.reserve(n);

    for (int i = 0; i < n; i++){
        if (population[i].dominatedCount == 0){
            current.push_back(i);
        }
    }

    while(!current.empty()){
        fronts_idx.push_back(current);
        vector<int> next;
        next.reserve(n);

        for (int p : current){
            for (int q : population[p].dominates){
                population[q].dominatedCount--;
                if (population[q].dominatedCount == 0){
                    next.push_back(q);
                }
            }
        }
        current = move(next);
    }

    Frontiers fronts;
    fronts.reserve(fronts_idx.size());
    int rank = 1;
    for (const auto& f : fronts_idx){
        Population layer;
        layer.reserve(f.size());
        for (int idx : f){
            population[idx].rank = rank;
            layer.push_back(population[idx]); 
        } 
        fronts.push_back(move(layer));
        rank++;
    }

    return fronts;
}

void sort_population(vector<Individual> &population){ // sort by dominance, then crowding distance
    sort(population.begin(), population.end(),
        [](const Individual& a, const Individual& b){
            // 1) Dominations has priority

            if (a.rank < b.rank) return true;   // a comes before b
            if (b.rank < a.rank) return false;  // b comes before a

            // 2) uses crowding distance as second
            if (!nearly_equal(a.crowdingDistance, b.crowdingDistance)) {
                return a.crowdingDistance > b.crowdingDistance;
            }

            // 3) make sure that if the previous way didn't sort, we still have determinism
            if (!nearly_equal(a.expectedReturn, b.expectedReturn)) {
                return a.expectedReturn > b.expectedReturn; // greatest return first
            }
            if (!nearly_equal(a.risk, b.risk)) {
                return a.risk < b.risk; // lower risk first
            }
            return a.picked.size() > b.picked.size(); // bigger portfolio size first
        }
    );
}

/**
 * Do operations to ensure the constraints
 */
void repair_population(vector<Individual> &population){ 

}



static int tournament_select_index(const vector<Individual>& pop, mt19937& rng) {
    const int n = (int)pop.size();
    int a = randi(rng, 0, n - 1);
    int b = randi(rng, 0, n - 1);
    while (b == a && n > 1) b = randi(rng, 0, n - 1);

    const Individual& A = pop[a];
    const Individual& B = pop[b];

    if (A.rank < B.rank) return a;
    if (B.rank < A.rank) return b;

    // higher crowding better
    if (A.crowdingDistance > B.crowdingDistance) return a;
    if (B.crowdingDistance > A.crowdingDistance) return b;

    // tie-break random
    return (rand01(rng) < 0.5) ? a : b;
}



void calculate_crowding_distance(Frontiers& frontiers) {
    const double INF = numeric_limits<double>::infinity();

    for (auto& frontier : frontiers) {
        const int m = (int)frontier.size();
        if (m == 0) continue;

        // Reset distances
        for (auto& ind : frontier) ind.crowdingDistance = 0.0;

        // If only 1 or 2 individuals, all are boundary points
        if (m <= 2) {
            for (auto& ind : frontier) ind.crowdingDistance = INF;
            continue;
        }

        // ----- Objective 1: expectedReturn (maximize) -----
        {
            vector<int> idx(m);
            for (int i = 0; i < m; i++) idx[i] = i;

            sort(idx.begin(), idx.end(),
                      [&](int a, int b) { return frontier[a].expectedReturn < frontier[b].expectedReturn; });

            const double fmin = frontier[idx.front()].expectedReturn;
            const double fmax = frontier[idx.back()].expectedReturn;

            // Boundary points
            frontier[idx.front()].crowdingDistance = INF;
            frontier[idx.back()].crowdingDistance  = INF;

            if (!nearly_equal(fmax, fmin)) {
                for (int k = 1; k < m - 1; k++) {
                    int i = idx[k];
                    if (isinf(frontier[i].crowdingDistance)) continue;

                    const double prev = frontier[idx[k - 1]].expectedReturn;
                    const double next = frontier[idx[k + 1]].expectedReturn;

                    frontier[i].crowdingDistance += (next - prev) / (fmax - fmin);
                }
            }
        }

        // ----- Objective 2: risk (minimize) -----
        // Crowding distance always uses sorted order (ascending).
        // Minimization vs maximization doesn't change the distance formula
        // (it depends on neighbor gaps), only the sort direction. Ascending is fine.
        {
            vector<int> idx(m);
            for (int i = 0; i < m; i++) idx[i] = i;

            sort(idx.begin(), idx.end(),
                      [&](int a, int b) { return frontier[a].risk < frontier[b].risk; });

            const double fmin = frontier[idx.front()].risk;
            const double fmax = frontier[idx.back()].risk;

            // Boundary points
            frontier[idx.front()].crowdingDistance = INF;
            frontier[idx.back()].crowdingDistance  = INF;

            if (!nearly_equal(fmax, fmin)) {
                for (int k = 1; k < m - 1; k++) {
                    int i = idx[k];
                    if (isinf(frontier[i].crowdingDistance)) continue;

                    const double prev = frontier[idx[k - 1]].risk;
                    const double next = frontier[idx[k + 1]].risk;

                    frontier[i].crowdingDistance += (next - prev) / (fmax - fmin);
                }
            }
        }
    }
}


int main(){
    PortfolioData data = PortfolioDataLoader::load_from_file("port1.txt");

    Population population = generate_population(POP_SIZE, data);
    population.reserve(POP_SIZE * 2);


    int generation = 0;

    Population nextGeneration;
    nextGeneration.reserve(population.size());

    //  if(DEBUG){
    //     cout << "Loaded " << data.mean.size() << " assets" << endl;
    //     cout << "Generate " << population.size() << " individuals" << endl;
    //     cout << "frontiers number " << frontiers.size() << endl;
    //     print_frontiers(frontiers, population);
    // }

    for(int i=0; i<GENERATIONS; i++){
        Frontiers frontiers = non_dominant_sort(population); 

        calculate_crowding_distance(frontiers);

        int frontierIndex=0;
        while(nextGeneration.size() + frontiers[frontierIndex].size() < POP_SIZE){
            nextGeneration.insert(nextGeneration.end(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].end());
            frontierIndex++;
        }
        sort_population(frontiers[frontierIndex]);
        int remainingSpace =  POP_SIZE-nextGeneration.size();
        nextGeneration.insert(nextGeneration.begin(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].begin() + remainingSpace);

        // Population offspring = generate_population(nextGeneration);
        // repair_population(offspring);

        // concatenation of nextGeneration + offspring
        // population.clear();
        // population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
        // population.insert(population.end(), offspring.begin(), offspring.end());
    }
   
}