#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <cstddef>

struct PortfolioData {
    int n = 0;
    std::vector<double> mean;   // size n
    std::vector<double> stdev;  // size n
    std::vector<double> corr;   // size n*n, row-major: corr[i*n + j]

    double corr_at(int i, int j) const {
        return corr[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)];
    }

    void set_corr(int i, int j, double v) {
        corr[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
             static_cast<std::size_t>(j)] = v;
    }

    // Covariância: cov_ij = sigma_i * sigma_j * rho_ij
    double cov_at(int i, int j) const {
        return stdev[static_cast<std::size_t>(i)] *
               stdev[static_cast<std::size_t>(j)] *
               corr_at(i, j);
    }
};

class PortfolioDataLoader {
public:
    // Lança std::runtime_error em caso de erro (mensagem bem explicativa)
    static PortfolioData load_from_file(const std::string& filepath);

private:
    static void validate_indices(int i1, int j1, int n, const std::string& filepath);
};
