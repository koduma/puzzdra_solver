# Puzzle & Dragons Solver

## Compiler:MinGW-W64

## Command 

- g++ -O2 -std=c++11 -fopenmp puzzdra_solver.cpp -o puzzdra_solver  
- g++ -O2 -std=c++11 -fopenmp -mbmi2 puzzdra_solver_BBver.cpp -o puzzdra_solver_BBver

# Benchmark

## Parameters

- Board = 5x6
- Drop = 6 (fire, water, wood, light, dark, heal)
- Direction = 4 (up,down,left,right)
- Problem = 10000
- Beam_width = 10000
- Max_turn = Trn = 150


| Version | CPU | Average Time | Average Combo |
| --- | --- | --- | --- |
| puzzdra_solver.cpp | Intel(R) Core(TM) i5-8250U | 0.2494s | 7.9907/7.9907 |
| puzzdra_solver_BBver.cpp | Intel(R) Core(TM) i5-8250U| 0.1805s | 7.9948/7.9948 |

<img src="https://user-images.githubusercontent.com/47982907/101321654-0b96e900-38a9-11eb-9c70-a8d9fa3d491d.jpg" width="300px" height="500px">
