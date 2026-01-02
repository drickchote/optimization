#include <iostream>
#include <list>
#include <tuple>

using namespace std;


// max return
// min risk

struct Individual{
    double expectedReturn;
    double risk;
    int rank;
    int dominatedCount;
    list<Individual> dominates;
};

list<Individual> generate_population(){

}

bool dominates(const Individual& individual1, const Individual& individual2){
    return (individual1.expectedReturn > individual2.expectedReturn &&  individual1.risk >= individual2.risk)
            || (individual1.expectedReturn >= individual2.expectedReturn &&  individual1.risk > individual2.risk);
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



int main(){
    list<Individual> population = generate_population();

    
    for(auto individual: population){
        if(individual.dominatedCount == 0){

        }
    }
}