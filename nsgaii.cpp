#include <iostream>
#include <list>
#include <tuple>

using namespace std;


static constexpr int POP_SIZE = 40;
static constexpr int GENERATIONS = 500;
static constexpr int NUMBER_OF_ASSETS = 31; // get from the file portx.txt
static constexpr double LB = 0.01; //
static constexpr double UB = 1.0; // 
static constexpr double K = 10; // max picked assets 



// max return
// min risk

struct Individual{
    vector<bool> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
    int rank;
    int dominatedCount;
    vector<Individual> dominates;
};

vector<double> generateDecision(int size){

}

vector<Individual> generate_population(int size){ // generates a random population
    vector<Individual> population;
    population.reserve(POP_SIZE);

    for(int i=0; i<size; i++){
        Individual individual = Individual();
        
    }
}

vector<Individual> generate_population(vector<Individual> population){ // generates a population with genetic operations

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
void repairPopulation(vector<Individual> &population){ 

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
        repairPopulation(offspring);

        // concatenation of nextGeneration + offspring
        population.clear();
        population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
        population.insert(population.end(), offspring.begin(), offspring.end());
    }
    
   
}