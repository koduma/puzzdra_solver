# Puzzle & Dragons Solver

## Compiler:MinGW-W64

## Compile Command 

- g++ -O2 -std=c++11 -fopenmp benchmark.cpp -o benchmark  
- g++ -O2 -std=c++11 -fopenmp -mbmi2 benchmark_BitBoard_ver.cpp -o benchmark_BitBoard_ver


# Benchmark

## Parameters

- Board = 5x6
- Drop = 6 (fire, water, wood, light, dark, heal)
- Direction = 4 (up,down,left,right)
- Problem = 1000
- Beam_width = 10000
- Max_turn = Trn = 150


| Version | CPU | Average Time | Average Combo |
| --- | --- | --- | --- |
| benchmark.cpp | Intel(R) Core(TM) i5-8250U | 0.295s | 7.991/7.991 |
| benchmark_BitBoard_ver.cpp | Intel(R) Core(TM) i5-8250U| 0.222s | 7.990/7.990 |

<img src="https://user-images.githubusercontent.com/47982907/101321654-0b96e900-38a9-11eb-9c70-a8d9fa3d491d.jpg" width="300px" height="500px">
