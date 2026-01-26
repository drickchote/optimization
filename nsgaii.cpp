#include <iostream>
#include <list>
#include <tuple>
#include <random>
#include "portfolio_data.hpp"

using namespace std;

static int DEBUG = 1;


static constexpr int POP_SIZE = 40;
static constexpr int GENERATIONS = 500;
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
using Frontiers = vector<vector<int>>; 

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
        for(auto individualIdx : frontier){
            print_individual(population.at(individualIdx));
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

    return std::sqrt(variance);
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

Population crossover(Individual individual1, Individual individual2){
    
}




bool dominates(const Individual& individual1, const Individual& individual2){
    return (individual1.expectedReturn > individual2.expectedReturn &&  individual1.risk <= individual2.risk)
            || (individual1.expectedReturn >= individual2.expectedReturn &&  individual1.risk < individual2.risk);
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

Frontiers non_dominant_sort(Population & population){
    evaluate(population);

    Frontiers fronts;
    vector<int> current;

    const int n = (int)population.size();
    current.reserve(n);

    for (int i = 0; i < n; i++){
        if (population[i].dominatedCount == 0){
            current.push_back(i);
        }
    }

    while(!current.empty()){
        fronts.push_back(current);
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

        current = std::move(next);
    }

    return fronts;
}

void sort(vector<Individual> &population){ // sort by crowding distance

}

/**
 * Do operations to ensure the constraints
 */
void repair_population(vector<Individual> &population){ 

}

static inline int randi(mt19937& rng, int lo, int hi_inclusive) {
    uniform_int_distribution<int> d(lo, hi_inclusive);
    return d(rng);
}

static inline double rand01(mt19937& rng) {
    uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng);
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




int main(){
    PortfolioData data = PortfolioDataLoader::load_from_file("port1.txt");

    Population population = generate_population(POP_SIZE, data);
    population.reserve(POP_SIZE * 2);


    int generation = 0;

    Population nextGeneration;
    nextGeneration.reserve(population.size());
    Frontiers frontiers = non_dominant_sort(population); 

     if(DEBUG){
        cout << "Loaded " << data.mean.size() << " assets" << endl;
        cout << "Generate " << population.size() << " individuals" << endl;
        cout << "frontiers number " << frontiers.size() << endl;
        print_frontiers(frontiers, population);
    }

    // for(int i=0; i<GENERATIONS; i++){
    //     vector<vector<Individual>> frontiers = non_dominant_sort(population); 

    //     int frontierIndex=0;
    //     while(nextGeneration.size() + frontiers[frontierIndex].size() < POP_SIZE){
    //         nextGeneration.insert(nextGeneration.end(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].end());
    //         frontierIndex++;
    //     }
    //     sort(frontiers[frontierIndex]);
    //     int remainingSpace =  POP_SIZE-nextGeneration.size();
    //     nextGeneration.insert(nextGeneration.begin(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].begin() + remainingSpace);

    //     vector<Individual> offspring = generate_population(nextGeneration);
    //     repair_population(offspring);

    //     // concatenation of nextGeneration + offspring
    //     population.clear();
    //     population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
    //     population.insert(population.end(), offspring.begin(), offspring.end());
    // }
    
   
}