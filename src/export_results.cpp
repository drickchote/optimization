#include "export_results.hpp"

#include <filesystem>
#include <cctype>
#include <stdexcept>

std::string portfolio_label_from_file(const std::string& portfolio_file) {
    const std::filesystem::path path(portfolio_file);
    std::string name = path.stem().string();

    if (name.size() >= 4 && name.rfind("port", 0) == 0) {
        return "Port " + name.substr(4);
    }

    if (!name.empty()) {
        name[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(name[0])));
    }

    return name;
}

std::string build_results_filepath(
    const std::string& algorithm,
    const std::string& portfolio_file,
    int k
) {
    const std::filesystem::path portfolio_path(portfolio_file);
    const std::string portfolio_name = portfolio_path.stem().string();

    std::string algorithm_dir = algorithm;
    for (char& character : algorithm_dir) {
        character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }

    const std::filesystem::path output_path =
        std::filesystem::path("runs") / algorithm_dir / (portfolio_name + "_k" + std::to_string(k) + ".csv");

    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    return output_path.string();
}

std::string build_results_header(
    const std::string& portfolio_file,
    int k,
    const std::string& algorithm
) {
    return portfolio_label_from_file(portfolio_file)
        + " - K = "
        + std::to_string(k)
        + " - "
        + algorithm;
}
