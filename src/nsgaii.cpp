#include <iostream>
#include <list>
#include <tuple>
#include <random>
#include <algorithm>
#include "portfolio_data.hpp"
#include "nsgaii.hpp"
#include <vector>
#include <iomanip>

using namespace std;

static int OUTPUT = 2; // 0 - None | 1 - Last Frontier: Risk Return | 2 - Each: Generation,Risk,Return

static constexpr int POP_SIZE = 250;
static constexpr int GENERATIONS = 400;
static constexpr double LB = 0.01; //
static constexpr double UB = 1.0; // 
static constexpr double K = 10; // max picked assets 
int NUMBER_OF_ASSETS = 0;

mt19937 rng(2);

using Population = NSGAII_Population;
using Frontiers = NSGAII_Frontiers; 
using Individual = NSGAII_Individual;
using Archive = vector<Individual>;


static inline int randi(mt19937& rng, int lo, int hi_inclusive) {
    uniform_int_distribution<int> d(lo, hi_inclusive);
    return d(rng);
}

static inline double rand01(mt19937& rng) {
    uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng);
}

static inline double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

void print_individual(Individual individual){
    if(OUTPUT){
        cout << setprecision(17) << individual.risk << " " << individual.expectedReturn;
    }

    cout << endl;
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

void print_pareto(Frontiers frontiers, Population population){
    int frontierNumber = 1;
    auto frontier = frontiers[0];
    frontierNumber++;
    for(auto individual : frontier){
        print_individual(individual);
    }
    cout << "------------" << endl;
}


double calculate_expected_return(const Individual& ind, const PortfolioData& data) {
    const int N = data.n;
    double er = 0.0;

    for (int i = 0; i < N; ++i) {
        er += ind.picked[i] * ind.weights[i] * data.mean[i];
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

    return variance;
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

/** Crossover operation based on (streichert, 2004b) */
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

    // 2) weights:
    for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
        const double p1 = individual1.weights[i];
        const double p2 = individual2.weights[i];

        child1.weights[i] = blx_gene(p1, p2);
        child2.weights[i] = blx_gene(p1, p2);
    }

    children.push_back(move(child1));
    children.push_back(move(child2));
    return children;
}

/** Mutation based on (streichert, 2004b)  */
void mutate(Individual& individual, double pm_bits = 0.1, double pm_real = 1.0, double sigma = 0.05) {
    // --- Bit-string one-point mutation (flip 1 gene) ---
    if (!individual.picked.empty() && rand01(rng) <= pm_bits) {
        int idx = randi(rng, 0, (int)individual.picked.size() - 1);

        // picked should be 0/1
        individual.picked[idx] = (individual.picked[idx] == 0) ? 1 : 0;
    }

    // --- Real-valued gaussian mutation ---
    // pm_real = 1.0 in the paper => always mutates (so this "if" always passes)
    if (!individual.weights.empty() && rand01(rng) <= pm_real) {
        normal_distribution<double> gauss(0.0, sigma);

        for (double& w : individual.weights) {
            w = clamp01(w + gauss(rng)); // keep genotype bounded (repair/normalization is separate)
        }
    }
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
    }

    // Generate weights for the selected assets
    for (int i = 0; i < size; ++i) {
        if (individual.picked[i] == 1) {
            individual.weights[i] = realDist(rng);
        }
    }

    return individual;
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
 * Repair only for constraint (Streichert et al.):
 *  - Keep only the K largest weights among selected assets
 *  - Set the rest to picked=0 and weight=0
 *  - Normalize remaining picked weights to sum to 1
 *  - Ensure at least one asset is selected
 *
 * This version is Lamarckian (it changes genotype: picked + weights).
 */
void repair_individual(Individual &individual) {
    auto ensure_at_least_one_picked = [&](Individual& ind) {
        int pickedCount = 0;
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) pickedCount += (ind.picked[i] != 0); 
        if (pickedCount > 0) return;
        // pick the max-weight gene (or random) and force it
        int best = 0;
        double bestW = -1.0;
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (ind.weights[i] > bestW) { bestW = ind.weights[i]; best = i; }
        }
        ind.picked[best] = 1;
        ind.weights[best] = 1.0;
    };

    auto apply_cardinality_keep_k_largest = [&](Individual& ind) {
        // First enforce consistency: if not picked -> weight 0; if picked -> clamp to [0,1]
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            ind.picked[i] = (ind.picked[i] != 0) ? 1 : 0;
            if (ind.picked[i] == 0) ind.weights[i] = 0.0;
            else ind.weights[i] = clamp01(ind.weights[i]);
        }

        // Collect indices of picked assets
        vector<int> idx;
        idx.reserve(NUMBER_OF_ASSETS);
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (ind.picked[i] == 1) idx.push_back(i);
        }

        if ((int)idx.size() <= (int)K) return;

        // Keep only top-K by weight among picked
        nth_element(
            idx.begin(),
            idx.begin() + (int)K,
            idx.end(),
            [&](int a, int b) { return ind.weights[a] > ind.weights[b]; }
        );

        vector<char> keep(NUMBER_OF_ASSETS, 0);
        for (int t = 0; t < (int)K; ++t) keep[idx[t]] = 1;

        // Zero out everything not kept
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (ind.picked[i] == 1 && !keep[i]) {
                ind.picked[i] = 0;
                ind.weights[i] = 0.0;
            }
        }
    };

    auto normalize = [&](Individual& ind) {
        double sum = 0.0;
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (ind.picked[i] == 0){
                ind.weights[i] = 0.0;
            } else {
                ind.weights[i] = clamp01(ind.weights[i]);
                sum += ind.weights[i];
            }

        }

        if (nearly_equal(sum, 0.0)) {
            // fallback: all mass on first picked asset
            int idx = -1;
            for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
                if (ind.picked[i]) { idx = i; break; }
            }
            for (int i = 0; i < NUMBER_OF_ASSETS; ++i) ind.weights[i] = 0.0;
            if (idx >= 0) ind.weights[idx] = 1.0;
            return;
        }

        int pickedCount = 0;
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) pickedCount += (ind.picked[i] != 0); 

        /* Adding lower bound constraint
        * Anagnostopoulos, 2011, pg. 4, section 3.2.2
        */
        for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
            if (ind.picked[i] == 1) {
                if(individual.weights[i] < 0){
                    cout << "Deu ruim vei"  << endl;
                }
                ind.weights[i] = LB  + (ind.weights[i] / sum) * ( 1 - LB * pickedCount);
            }
        }
    };



    // // 1) ensure at least one selected
    ensure_at_least_one_picked(individual);

    // 2) cardinality repair (keep only K largest)
    apply_cardinality_keep_k_largest(individual);

    // 3) normalize to sum 1
    normalize(individual);

    // ) final consistency
    for (int i = 0; i < NUMBER_OF_ASSETS; ++i) {
        if (!individual.picked[i]) individual.weights[i] = 0.0;
        
    }
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



bool dominated_or_existing_solution(Archive& archive, const Individual& ind){
    for (const auto& ub : archive){
        if (dominates(ub, ind) || (nearly_equal(ub.expectedReturn, ind.expectedReturn) && nearly_equal(ub.risk, ind.risk))){
            return true;
        }
    }
    return false;
}


void add_to_archive(Archive& archive, Individual& individual){
    if(dominated_or_existing_solution(archive, individual)){
        return;
    }

    for (auto it = archive.begin(); it != archive.end(); ) {
        if (dominates(individual, *it)) {
            it = archive.erase(it);  
        } else {
            ++it;
        }
    }
    archive.push_back(individual);
}

void add_population_to_archive(Archive& archive, Population& population){
    for(auto individual : population){
        add_to_archive(archive, individual);
    }
}

/* Generates a random population*/
Population generate_population(int size, PortfolioData data){
    Population population;
    population.reserve(POP_SIZE);

    for(int i=0; i<size; i++){
        Individual individual = generate_decision(NUMBER_OF_ASSETS, rng);
        repair_individual(individual);
        individual.expectedReturn = calculate_expected_return(individual, data);
        individual.risk = calculate_risk(individual, data);
        population.push_back(individual);
    }

    return population;
}

/* Generates a population with genetic operations */
Population generate_population(Population population, PortfolioData data){
    Population offspring;
    offspring.reserve(POP_SIZE);

    while ((int)offspring.size() < POP_SIZE) {
        int p1 = tournament_select_index(population, rng);
        int p2 = tournament_select_index(population, rng);
        while (p2 == p1 && (int)population.size() > 1) {
            p2 = tournament_select_index(population, rng);
        }

        Population children = crossover(population[p1], population[p2]);

        for (auto& child : children) {
            // mutation + repair/normalize
            mutate(child);
            repair_individual(child);
            child.risk = calculate_risk(child, data);
            child.expectedReturn = calculate_expected_return(child, data);

            offspring.push_back(move(child));
            if ((int)offspring.size() >= POP_SIZE) break;
        }
    }

    return offspring;
}

void print_csv_population(Archive archive, int generation){
    for(auto i : archive){
        cout << generation << "," << setprecision(17) << i.risk << "," << i.expectedReturn << endl;
    }
}



Population run_nsgaII(){
    Archive archive = {};
    archive.reserve(POP_SIZE * GENERATIONS);
    PortfolioData data = PortfolioDataLoader::load_from_file("port1.txt"); 
    NUMBER_OF_ASSETS = data.mean.size();
    portfolioData = data;
    Population population = generate_population(POP_SIZE, data);
    Frontiers frontiers = non_dominant_sort(population); 
    
    population.reserve(POP_SIZE * 2);
    Population nextGeneration;
    nextGeneration.reserve(POP_SIZE*2);
    
    if(OUTPUT == 2){
        cout << "generation,risk,expectedReturn" << endl;
    }

    for(int i=0; i<GENERATIONS; i++){
        nextGeneration.clear();
        // Torneio binário: Selecionar 2 indivíduos, eles competem no primeiro objetivo: Eles competem no primeiro objetivo
        // Seleciona 2 outros indivíduos e eles competem no segundo objetivo
        
        calculate_crowding_distance(frontiers);

        int frontierIndex=0;
        while(nextGeneration.size() + frontiers[frontierIndex].size() < POP_SIZE){
            nextGeneration.insert(nextGeneration.end(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].end());
            frontierIndex++;
        }
        sort_population(frontiers[frontierIndex]);
        int remainingSpace =  POP_SIZE-nextGeneration.size();
        nextGeneration.insert(nextGeneration.end(), frontiers[frontierIndex].begin(), frontiers[frontierIndex].begin() + remainingSpace);

        Population offspring = generate_population(nextGeneration, data);
        
        // concatenation of nextGeneration + offspring
        population.clear();
        population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
        population.insert(population.end(), offspring.begin(), offspring.end());

        frontiers = non_dominant_sort(population); 
        // add_population_to_archive(archive, frontiers[0]);

        if(OUTPUT == 2){
            print_csv_population(frontiers[0], i+1);
        }
    }

    if (OUTPUT == 1){
        print_population(frontiers[0]);
    } 

    add_population_to_archive(archive, frontiers[0]);
    return archive;
}


#ifndef BB
int main(){
    Population archive = run_nsgaII();
    return 0;
}
#endif