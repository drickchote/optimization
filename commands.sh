# Single-run builds (from src/)
g++ -std=c++17 nsgaii.cpp portfolio_data.cpp export_results.cpp experiment_config.cpp experiment_cli.cpp -o nsgaii && ./nsgaii

# For Apple Silicon Macs
g++ -std=c++17 nsgaii.cpp portfolio_data.cpp export_results.cpp experiment_config.cpp experiment_cli.cpp \
  -I/opt/homebrew/include \
  -L/opt/homebrew/lib \
  -o nsgaii -lpagmo && ./nsgaii --portfolio ../inputs/port1.txt --k 10


g++ -std=c++17 branch_and_bound.cpp portfolio_data.cpp nsgaii.cpp export_results.cpp experiment_config.cpp experiment_cli.cpp -o branch_and_bound && ./branch_and_bound


g++ -std=c++17 -DBB branch_and_bound.cpp portfolio_data.cpp nsgaii.cpp export_results.cpp experiment_config.cpp experiment_cli.cpp \
  -I/opt/homebrew/include \
  -L/opt/homebrew/lib \
  -I/Library/gurobi1203/macos_universal2/include \
  -L/Library/gurobi1203/macos_universal2/lib \
  -lgurobi_c++ -lgurobi120 \
  -lm -lpthread \
  -o branch_and_bound -lpagmo && ./branch_and_bound --portfolio ../inputs/port4.txt --k 10 - -


g++ -std=c++17 branch_and_bound.cpp portfolio_data.cpp nsgaii.cpp export_results.cpp experiment_config.cpp experiment_cli.cpp -o branch_and_bound \
  -I/opt/homebrew/include \
  -L/opt/homebrew/lib \
  -o branch_and_bound -lpagmo && ./branch_and_bound --portfolio ../inputs/port1.txt --k 10 - -


# Batch sweep (all portfolios and K values from experiments.json)
cd src && python3 scripts/run_experiments.py

g++ normalize.cc -o normalize && ./normalize normalize_param.txt boundfile.txt ~/learning/mestrado/statistical_tests_suit/algorithm_results/NSGA2/port2/run_1.txt normalized_nsga2_port2
