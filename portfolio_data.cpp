#include "portfolio_data.hpp"

#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

static std::runtime_error make_err(const std::string& filepath, int lineNo, const std::string& msg) {
    std::ostringstream oss;
    if (lineNo > 0) oss << filepath << ":" << lineNo << ": ";
    else oss << filepath << ": ";
    oss << msg;
    return std::runtime_error(oss.str());
}

void PortfolioDataLoader::validate_indices(int i1, int j1, int n, const std::string& filepath) {
    if (i1 < 1 || i1 > n || j1 < 1 || j1 > n) {
        std::ostringstream oss;
        oss << "Correlation indices out of range: (" << i1 << ", " << j1 << "). Expected 1..N (N=" << n << ").";
        throw std::runtime_error(filepath + ": " + oss.str());
    }
}

PortfolioData PortfolioDataLoader::load_from_file(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }

    PortfolioData d;

    std::string line;
    int lineNo = 0;

    auto next_non_empty_line = [&]() -> bool {
        while (std::getline(in, line)) {
            ++lineNo;
            bool allSpace = true;
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) { allSpace = false; break; }
            }
            if (!allSpace) return true;
        }
        return false;
    };

    if (!next_non_empty_line()) {
        throw make_err(filepath, lineNo, "File is empty; expected N on the first line.");
    }
    {
        std::istringstream iss(line);
        if (!(iss >> d.n) || d.n <= 0) {
            throw make_err(filepath, lineNo, "Failed to read N (number of assets) or N <= 0.");
        }
    }

    const std::size_t N = static_cast<std::size_t>(d.n);
    d.mean.assign(N, 0.0);
    d.stdev.assign(N, 0.0);
    d.corr.assign(N * N, 0.0);

    for (std::size_t i = 0; i < N; ++i) {
        d.corr[i * N + i] = 1.0;
    }

    for (int i = 0; i < d.n; ++i) {
        if (!next_non_empty_line()) {
            throw make_err(filepath, lineNo, "Unexpected EOF while reading mean/stdev section.");
        }

        double mu = 0.0, sd = 0.0;
        std::istringstream iss(line);
        if (!(iss >> mu >> sd)) {
            throw make_err(filepath, lineNo, "Malformed mean/stdev line. Expected: <mean> <stdev>.");
        }
        if (sd < 0.0) {
            throw make_err(filepath, lineNo, "Standard deviation cannot be negative.");
        }

        d.mean[static_cast<std::size_t>(i)] = mu;
        d.stdev[static_cast<std::size_t>(i)] = sd;
    }

    while (next_non_empty_line()) {
        int i1 = 0, j1 = 0;
        double rho = 0.0;

        std::istringstream iss(line);
        if (!(iss >> i1 >> j1 >> rho)) {
            throw make_err(filepath, lineNo, "Malformed correlation line. Expected: <i> <j> <corr>.");
        }

        validate_indices(i1, j1, d.n, filepath);

        if (rho < -1.0000001 || rho > 1.0000001) {
            throw make_err(filepath, lineNo, "Correlation value seems out of range [-1, 1].");
        }

        const int i = i1 - 1;
        const int j = j1 - 1;

        d.set_corr(i, j, rho);
        d.set_corr(j, i, rho); // simetria
    }

    return d;
}
