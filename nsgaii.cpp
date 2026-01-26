#include <iostream>
#include <list>
#include <tuple>
#include <random>

using namespace std;


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
    int dominatedCount;
    vector<Individual> dominates;
    double crowdingDistance = 0.0;
};

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

    // at least one should be selected
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

vector<Individual> generate_population(int size){ // generates a random population
    vector<Individual> population;
    population.reserve(POP_SIZE);

    for(int i=0; i<size; i++){
        Individual individual = generate_decision(NUMBER_OF_ASSETS, rng);
        
    }
}

vector<Individual> generate_population(vector<Individual> population){ // generates a population with genetic operations

}

vector<Individual> crossover(Individual individual1, Individual individual2){
    
}




bool dominates(const Individual& individual1, const Individual& individual2){
    return (individual1.expectedReturn > individual2.expectedReturn &&  individual1.risk <= individual2.risk)
            || (individual1.expectedReturn >= individual2.expectedReturn &&  individual1.risk < individual2.risk);
}


void evaluate(list<Individual> &population){
    for(auto individual1 : population){
        for(auto individual2 : population){
            if(dominates(individual1, individual2)){
                individual2.dominatedCount++;
                individual1.dominates.push_back(individual2);
            } else if(dominates(individual2, individual1)){
                individual1.dominatedCount++;
                individual2.dominates.push_back(individual1);
            }
        }
    }
}

vector<vector<Individual>> non_dominant_sort(vector<Individual> population){

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
    vector<Individual> population = generate_population(POP_SIZE);
    population.reserve(POP_SIZE * 2);

    int generation = 0;
    
    vector<Individual> nextGeneration;
    nextGeneration.reserve(population.size());
    
    
    for(int i=0; i<GENERATIONS; i++){
        vector<vector<Individual>> frontiers = non_dominant_sort(population); 

        int frontierIndex=0;
        while(nextGeneration.size() + frontiers[frontierIndex].size() < POP_SIZE){
            nextGeneration.insert(nextGeneration.end(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].end());
            frontierIndex++;
        }
        sort(frontiers[frontierIndex]);
        int remainingSpace =  POP_SIZE-nextGeneration.size();
        nextGeneration.insert(nextGeneration.begin(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].begin() + remainingSpace);

        vector<Individual> offspring = generate_population(nextGeneration);
        repair_population(offspring);

        // concatenation of nextGeneration + offspring
        population.clear();
        population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
        population.insert(population.end(), offspring.begin(), offspring.end());
    }
    
   
}