#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>

struct Individual {
    std::vector<int> picked;
    std::vector<double> weights;
    double expectedReturn = 0.0;
    double risk = 0.0;

    double get_obj(int objective) const {
        return objective == 0 ? risk : expectedReturn;
    }
};

#endif
