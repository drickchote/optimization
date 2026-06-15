#ifndef EXPERIMENT_CLI_HPP
#define EXPERIMENT_CLI_HPP

#include <string>

struct ExperimentOptions {
    std::string portfolio_file;
    int k = 0;
    std::string checkpoint_in;
    std::string checkpoint_out;
};

ExperimentOptions parse_experiment_cli(int argc, char** argv, bool allow_checkpoints);

void apply_experiment_options(const ExperimentOptions& options);

void print_nsgaii_usage(const char* program);

void print_branch_and_bound_usage(const char* program);

#endif
