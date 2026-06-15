#ifndef GRID_HPP
#define GRID_HPP

#include "param.h"
#include <cstring>

class Grid {
private:
    static constexpr int GRID_SIZE =
        1 << (NUMBER_OF_OBJECTIVES * GRID_DEPTH);

    int positions[GRID_SIZE];

public:
    Grid() {
        clear_grid();
    }

    int get_position_count(int position) const {
        if (position < 0 || position >= GRID_SIZE) {
            return -1;
        }

        return positions[position];
    }

    void add_grid(int position) {
        ++positions[position];
    }

    void remove_grid(int position) {
        --positions[position];
    }

    void clear_grid() {
        memset(positions, 0, sizeof(positions));
    }

    int get_size() const {
        return GRID_SIZE;
    }
};

#endif