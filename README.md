[![ubuntu-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/ubuntu-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/ubuntu-latest.yml)
[![windows-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/windows-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/windows-latest.yml)
[![macos-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/macos-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/macos-latest.yml)

# Contact https://twitter.com/tekitouk

# Puzzle & Dragons Solver

## Movie
[![](https://img.youtube.com/vi/eL6Sf42STzA/0.jpg)](https://www.youtube.com/watch?v=eL6Sf42STzA)

## Compiler:g++

## Compile Command

| OS | Version | Compile Command |
| --- | --- | --- |
|  Any | benchmark.cpp | g++ -O2 -std=c++11 -fopenmp benchmark.cpp -o benchmark  |
|  Any | benchmark_BitBoard_ver.cpp | g++ -O2 -std=c++11 -fopenmp -mbmi2 benchmark_BitBoard_ver.cpp -o benchmark_BitBoard_ver |
|  Linux | Excalibur.cpp | g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large -ldl  |
|  Windows | Excalibur.cpp | g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large  |
|  MacOS | Excalibur.cpp | g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -ldl  |

# Benchmark

## Parameters

- Board = 5x6
- Drop = 6 (fire, water, wood, light, dark, heal)
- Direction = 4 (up,down,right,left)
- Problem = 10000
- Beam_width = 10000
- Max_turn = Trn = 150


| Version | CPU | Average Time | Average Combo |
| --- | --- | --- | --- |
| benchmark.cpp | Intel(R) Core(TM) i5-8250U | 0.2301s | 8.0136/8.0136 |
| benchmark_BitBoard_ver.cpp | Intel(R) Core(TM) i5-8250U| 0.1618s | 7.9968/7.9968 |

<img src="https://user-images.githubusercontent.com/47982907/101321654-0b96e900-38a9-11eb-9c70-a8d9fa3d491d.jpg" width="300px" height="500px">
