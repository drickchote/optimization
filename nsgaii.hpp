#include <list>

using namespace std;

struct NSGAII_Individual{
    vector<int> picked; // list of assets weights
    vector<double> weights; // list of assets weights
    double expectedReturn; // maximize 
    double risk; // minimize
    int rank;
    int dominatedCount = 0;
    vector<int> dominates;
    double crowdingDistance = 0.0;
};

using NSGAII_Population = vector<NSGAII_Individual>;
using NSGAII_Frontiers = vector<vector<NSGAII_Individual>>; 

inline bool nearly_equal(double a, double b, double eps = 1e-12) {
    return fabs(a - b) <= eps;
}

NSGAII_Population run_nsgaII();

extern int NUMBER_OF_ASSETS;