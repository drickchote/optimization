#ifndef BOUNDED_PARETO_SET_CPP
#define BOUNDED_PARETO_SET_CPP

#include "pareto_set.cpp"

class BoundedParetoSet : public ParetoSet {
public:
    static constexpr int MAX_ARCHIVE_SIZE = 5000;
};

#endif
