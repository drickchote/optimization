#include <iostream>
#include <list>
#include <tuple>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"
#include <vector>

using namespace std;

using Individual = NSGAII_Individual;
using Population = NSGAII_Population;
using Archive = vector<Individual>;

const int OUTPUT = 1;




bool dominates(const Individual& a, const Individual& b){
    const bool betterReturn =
        (a.expectedReturn > b.expectedReturn);

    const bool betterRisk =
        (a.risk < b.risk);

    const bool noWorseReturn =
        (a.expectedReturn > b.expectedReturn) ||
        nearly_equal(a.expectedReturn, b.expectedReturn);

    const bool noWorseRisk =
        (a.risk < b.risk) ||
        nearly_equal(a.risk, b.risk);

    return (betterReturn && noWorseRisk)
        || (betterRisk && noWorseReturn);
}


bool dominated_or_existing_solution(Archive& archive, const Individual& ind){
    for (const auto& ub : archive){
        if (dominates(ub, ind) || (nearly_equal(ub.expectedReturn, ind.expectedReturn) && nearly_equal(ub.risk, ind.risk))){
            return true;
        }
    }
    return false;
}


void add_to_archive(Archive& arquive, Individual individual){
    if(dominated_or_existing_solution(arquive, individual)){
        return;
    }
    for (auto it = arquive.begin(); it != arquive.end(); ) {
        if (dominates(individual, *it)) {
            it = arquive.erase(it);  
        } else {
            ++it;
        }
    }

    arquive.push_back(individual);
}

void add_population_to_archive(Archive& arquive, Population& population){
    for(auto individual : population){
        add_to_archive(arquive, individual);
    }
}

void print_individual(Individual individual){
    if(OUTPUT){
        cout << individual.risk << " " << individual.expectedReturn << endl;
        return;
    }
    cout << "______________" << endl;
    cout << "Expected Return: " << individual.expectedReturn << endl;
    cout << "Risk: " << individual.risk << endl;
    cout << "______________" << endl;
}


void print_population(Population &population){
    for(auto i : population){
        print_individual(i);
    }
}

Archive make_population(vector<vector<double>> points){
    Archive archive = {};
    for(auto i : points){
        Individual individual = {};
        individual.risk = i[0];
        individual.expectedReturn = i[1];
        archive.push_back(individual);
    }

    return archive;
}

int main(){
    Archive arquive = {};
    Population population = make_population({
        {0.0344674,0.00389088},
        {0.0344678,0.00389088}
    });

    add_population_to_archive(arquive, population);
    print_population(arquive);
    return 0;
}