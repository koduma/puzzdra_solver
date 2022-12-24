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

//Excalibur.cppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/Excalibur.cpp

//hash_map.hpp,loguru.cpp,loguru.hppをダウンロード
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/hash_map.hpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.cpp
wget --no-check-certificate https://raw.githubusercontent.com/koduma/puzzdra_solver/master/loguru.hpp

//ビーム幅調整
vi Excalibur.cpp

//コンパイル
Linux:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large -ldl
Windows10,Windows11:g++ -O2 -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -mcmodel=large
MacOS:g++ -std=c++11 -fopenmp -mbmi2 -lpthread Excalibur.cpp loguru.cpp -o Excalibur -ldl

//run
./Excalibur

//input
*/
