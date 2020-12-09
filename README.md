# Puzzle & Dragons Solver

## Compiler:MinGW-W64

## Command 

- g++ -O2 -std=c++11 -fopenmp puzzdra_solver.cpp -o puzzdra_solver  
- g++ -O2 -std=c++11 -fopenmp -mbmi2 puzzdra_solver_BBver.cpp -o puzzdra_solver_BBver

# Benchmark

## Parameters

- 5x6 BORAD
- PROBLEM = 10000
- BEAM_WIDTH = 10000
- MAX_TURN = 150


| Version | CPU | Average Time | Average Combo |
| --- | --- | --- | --- |
| puzzdra_solver.cpp | Intel(R) Core(TM) i7-3770 | 0.2843s | 8.0010/8.0010 |
| puzzdra_solver_BBver.cpp | Intel(R) Core(TM) i5-8250U| 0.2130s | 8.0143/8.0143 |

<img src="https://user-images.githubusercontent.com/47982907/101321654-0b96e900-38a9-11eb-9c70-a8d9fa3d491d.jpg" width="300px" height="500px">
