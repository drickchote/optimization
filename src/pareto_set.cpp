#ifndef PARETO_SET_CPP
#define PARETO_SET_CPP

#include <vector>
#include <cmath>
#include <cassert>

#include "individual.hpp"
#include "param.h"
#include "grid.cpp"

struct Range {
    double min;
    double max;
};

class ParetoSet {
protected:
    Range new_range[NUMBER_OF_OBJECTIVES];
    Range current_range[NUMBER_OF_OBJECTIVES];
    Grid grid;

    int calculate_grid_position(const Individual& individual) const {
        int bit = 0;
        int grid_position = 0;

        for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
            double start = current_range[objective].min;
            double end = current_range[objective].max;
            double middle = (start + end) / 2.0;

            for (int depth = 0; depth < GRID_DEPTH; ++depth) {
                if (individual.get_obj(objective) >= middle) {
                    grid_position |= (1 << bit);
                    start = middle;
                } else {
                    end = middle;
                }

                middle = (start + end) / 2.0;
                ++bit;
            }
        }

        return grid_position;
    }

    void reset_ranges() {
        constexpr double INF = 1e9;

        for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
            current_range[objective].min = new_range[objective].min = INF;
            current_range[objective].max = new_range[objective].max = -INF;
        }
    }

    void update_grid(const std::vector<Individual>& archive) {
        grid.clear_grid();
        reset_ranges();

        for (const Individual& individual : archive) {
            for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
                current_range[objective].min =
                    new_range[objective].min =
                        std::min(current_range[objective].min, individual.get_obj(objective));

                current_range[objective].max =
                    new_range[objective].max =
                        std::max(current_range[objective].max, individual.get_obj(objective));
            }
        }

        for (const Individual& individual : archive) {
            grid.add_grid(calculate_grid_position(individual));
        }
    }

public:
    ParetoSet() {
        reset_ranges();
    }

    int get_position_count(const Individual& individual) const {
        return grid.get_position_count(calculate_grid_position(individual));
    }

    void rebuild(const std::vector<Individual>& archive) {
        update_grid(archive);
    }

    void remove_individual(const Individual& individual) {
        grid.remove_grid(calculate_grid_position(individual));
    }

    void add_individual(const Individual& individual) {
        for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
            if (new_range[objective].min > new_range[objective].max) {
                current_range[objective].min =
                    new_range[objective].min =
                        individual.get_obj(objective);

                current_range[objective].max =
                    new_range[objective].max =
                        individual.get_obj(objective);
            }
        }

        grid.add_grid(calculate_grid_position(individual));

        for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
            new_range[objective].min =
                std::min(new_range[objective].min, individual.get_obj(objective));

            new_range[objective].max =
                std::max(new_range[objective].max, individual.get_obj(objective));
        }
    }

    void finalize_addition(const std::vector<Individual>& archive) {
        for (int objective = 0; objective < NUMBER_OF_OBJECTIVES; ++objective) {
            const bool min_changed =
                std::fabs(new_range[objective].min - current_range[objective].min)
                > 0.1 * current_range[objective].min;

            const bool max_changed =
                std::fabs(new_range[objective].max - current_range[objective].max)
                > 0.1 * current_range[objective].max;

            if (min_changed || max_changed) {
                update_grid(archive);
                return;
            }
        }
    }

    bool check_grid(const std::vector<Individual>& archive) const {
        unsigned total = 0;

        for (int position = 0; position < grid.get_size(); ++position) {
            total += grid.get_position_count(position);
        }

        return total == archive.size();
    }
};

#endif
