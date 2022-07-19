[![ubuntu-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/ubuntu-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/ubuntu-latest.yml)
[![windows-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/windows-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/windows-latest.yml)
[![macos-latest](https://github.com/koduma/puzzdra_solver/actions/workflows/macos-latest.yml/badge.svg?branch=master)](https://github.com/koduma/puzzdra_solver/actions/workflows/macos-latest.yml)

# コンタクト https://twitter.com/tekitouk

# 懸賞問題

  手数を小さく出来た方には、1問につきitunesカードもしくはgoogle playカード1万円分差し上げます。

	layout=367254402726710107527213362754
	:手数=52,10combo
	layout=047631151072370164261053045210
	:手数=50,10combo
	layout=242242100331023100110324132543
	:手数=26,9combo
	layout=201053210251533425501353123221
	:手数=26,9combo
	layout=015315151020442313510540210411
	:手数=27,9combo
	layout=432015152244350331552132312515
	:手数=31,9combo
	layout=323243441332042002331313014300
	:手数=19,8combo
	layout=225530333313140355004550251403
	:手数=24,9combo
	layout=224234425402054400304510125043
	:手数=30,8combo
	layout=053241405407470557104053134522
	:手数=41,10combo

# Puzzle & Dragons Solver

## Movie

https://user-images.githubusercontent.com/47982907/147241543-7e7523d8-798a-4771-9f16-2c7dd535ac06.mp4



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
| benchmark.cpp | Intel(R) Core(TM) i5-8250U | 0.2501s | 8.0089/8.0089 |
| benchmark_BitBoard_ver.cpp | Intel(R) Core(TM) i5-8250U| 0.1530s | 7.9964/7.9964 |

<img src="https://user-images.githubusercontent.com/47982907/101321654-0b96e900-38a9-11eb-9c70-a8d9fa3d491d.jpg" width="300px" height="500px">
