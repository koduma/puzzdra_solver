/*
Windows10,Windows11,Linux,MacOS

Linux導入手続き

//メモリ容量確認
free -h

//g++インストール
sudo apt install -y g++

//wgetインストール
sudo apt-get update
sudo apt-get install -y wget

//test3.cppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/test3.cpp

//hash_map.hpp,loguru.cpp,loguru.hppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/hash_map.hpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.cpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.hpp

//ビーム幅調整
vi test3.cpp

//コンパイル
Linux:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread test3.cpp loguru.cpp -o test3 -mcmodel=large -ldl
Windows10,Windows11:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread test3.cpp loguru.cpp -o test3 -mcmodel=large
MacOS:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread test3.cpp loguru.cpp -o test3 -ldl

//run
./test3

//input
*/
