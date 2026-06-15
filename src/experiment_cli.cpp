#include "experiment_cli.hpp"

#include "experiment_config.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

ExperimentOptions parse_experiment_cli(int argc, char** argv, bool allow_checkpoints) {
    ExperimentOptions options;
    options.portfolio_file = PORTFOLIO_FILE;
    options.k = K;

    std::vector<std::string> positional;

    for (int index = 1; index < argc; ++index) {
        const char* arg = argv[index];

        if (std::strcmp(arg, "--portfolio") == 0) {
            if (index + 1 >= argc) {
                throw std::runtime_error("Missing value for --portfolio");
            }
            options.portfolio_file = argv[++index];
            continue;
        }

        if (std::strcmp(arg, "--k") == 0) {
            if (index + 1 >= argc) {
                throw std::runtime_error("Missing value for --k");
            }
            options.k = std::stoi(argv[++index]);
            continue;
        }

        if (arg[0] == '-' && arg[1] != '\0') {
            throw std::runtime_error(std::string("Unknown option: ") + arg);
        }

        positional.emplace_back(arg);
    }

    if (allow_checkpoints) {
        if (positional.size() == 2) {
            options.checkpoint_in = positional[0];
            options.checkpoint_out = positional[1];

            if (options.checkpoint_in == "-") {
                options.checkpoint_in.clear();
            }
            if (options.checkpoint_out == "-") {
                options.checkpoint_out.clear();
            }
        } else if (!positional.empty()) {
            throw std::runtime_error("Expected zero or two checkpoint positional arguments");
        }
    } else if (!positional.empty()) {
        throw std::runtime_error("Unexpected positional arguments");
    }

    return options;
}

void apply_experiment_options(const ExperimentOptions& options) {
    PORTFOLIO_FILE = options.portfolio_file;
    K = options.k;
}

void print_nsgaii_usage(const char* program) {
    std::cerr
        << "Usage: " << program << " [--portfolio path] [--k N]\n"
        << "  --portfolio  Portfolio input file (default: "
        << DEFAULT_PORTFOLIO_FILE << ")\n"
        << "  --k          Maximum number of assets in a portfolio (default: "
        << DEFAULT_K << ")\n";
}

void print_branch_and_bound_usage(const char* program) {
    std::cerr
        << "Usage: " << program << " [--portfolio path] [--k N] [checkpoint_in checkpoint_out]\n"
        << "  --portfolio     Portfolio input file (default: "
        << DEFAULT_PORTFOLIO_FILE << ")\n"
        << "  --k             Maximum number of assets in a portfolio (default: "
        << DEFAULT_K << ")\n"
        << "  checkpoint_in   File to resume from, or '-' to start from NSGA-II + B&B\n"
        << "  checkpoint_out  File written periodically, or '-' to disable\n";
}
