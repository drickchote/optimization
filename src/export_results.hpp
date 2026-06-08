#ifndef EXPORT_RESULTS_HPP
#define EXPORT_RESULTS_HPP

#include <fstream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

std::string portfolio_label_from_file(const std::string& portfolio_file);

std::string build_results_filepath(
    const std::string& algorithm,
    const std::string& portfolio_file,
    int k
);

std::string build_results_header(
    const std::string& portfolio_file,
    int k,
    const std::string& algorithm
);

template<typename Individual>
void export_results_csv(
    const std::string& filepath,
    const std::string& portfolio_file,
    int k,
    const std::string& algorithm,
    const std::vector<Individual>& individuals
) {
    const std::filesystem::path output_path(filepath);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    const std::string header = build_results_header(portfolio_file, k, algorithm);

    std::ofstream out(filepath);
    if (!out) {
        throw std::runtime_error("Failed to open results file: " + filepath);
    }

    out << header << '\n';
    out << std::setprecision(17);

    for (const Individual& individual : individuals) {
        out << individual.risk << ' ' << individual.expectedReturn << '\n';
    }
}

#endif
