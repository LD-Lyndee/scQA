g++ -c State.cpp
g++ -c Qualitative.cpp
g++ -c Quantitative.cpp
g++ -c Cluster.cpp
g++ -c main.cpp
g++ -o main main.o State.o Cluster.o Qualitative.o Quantitative.o
