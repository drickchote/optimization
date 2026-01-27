#include "portfolio_data.hpp"
#include <iostream>

int main() {
    auto data = PortfolioDataLoader::load_from_file("port1.txt");

    std::cout << "N=" << data.n << "\n";
    std::cout << "mean[0]=" << data.mean[0] << " stdev[0]=" << data.stdev[0] << "\n";
    std::cout << "corr(0,1)=" << data.corr_at(0,1) << "\n";
    std::cout << "cov(0,1)=" << data.cov_at(0,1) << "\n";
}