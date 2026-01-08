// Implementar SCH

#include <iostream>
#include <list>
#include <tuple>

using namespace std;


static constexpr int POP_SIZE = 40;
static constexpr int GENERATIONS = 500;
static constexpr int NUMBER_OF_ASSETS = 31; // get from the file portx.txt


// max return
// min risk

struct Individual{
    vector<double> decision; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
    int rank;
    int dominatedCount;
    vector<Individual> dominates;
};

vector<double> generateDecision(int size){

}

vector<Individual> generate_population(int size){
    vector<Individual> population;
    population.reserve(POP_SIZE);

    for(int i=0; i<size; i++){
        Individual individual = Individual();
        individual.decision = 
    }
}


bool dominates(const Individual& individual1, const Individual& individual2){
    return (individual1.expectedReturn > individual2.expectedReturn &&  individual1.risk <= individual2.risk)
            || (individual1.expectedReturn >= individual2.expectedReturn &&  individual1.risk < individual2.risk);
}

// Objective 1 - max expected
double expectedReturn(Individual individual){

}

// Objective 2
double risk(Individual individual){

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



int main(){
    vector<Individual> population = generate_population(POP_SIZE);
    
    int generation = 0;
    
    for(int i=0; i<GENERATIONS; i++){
        vector<vector<Individual>> frontiers = non_dominant_sort(population); 

    }
    
   
}