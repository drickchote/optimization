#ifndef PARAM_H
#define PARAM_H

// Knowles grid archive
#define NUMBER_OF_OBJECTIVES 2
#define GRID_DEPTH 5

// NSGA-II
#define POP_SIZE 250
#define GENERATIONS 400

// Portfolio problem
#define K 10
#define WEIGHT_LOWER_BOUND 0.01
#define WEIGHT_UPPER_BOUND 1.0

// Branch and bound
#define DEBUG 0

// Input data
#define PORTFOLIO_FILE "../inputs/port1.txt"

//#define ASSERT
#ifdef ASSERT
#define ASS(x) x
#else
#define ASS(x)
#endif

#endif
