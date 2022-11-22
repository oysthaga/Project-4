# Project-4
Project 4

MainCopy.cpp is for problems 4-6, while Main.cpp is parallelized. To compile Main.cpp:

g++ -o3  Main.cpp -larmadillo  -o Main.exe -fopenmp

export OMP_NUM_THREADS=8

to run Main.cpp:

./Main.exe [Number of Monte Carlo Cycles] [L] [0=unordered initial state, 1= ordered] [Number of temperature-steps] [minimum temperature] [maximum temperature]

To compile MainCopy.cpp:

g++ -o3  Main.cpp -larmadillo  -o Main.exe 

to run Main.cpp:

./Main.exe [Number of Monte Carlo Cycles] [temperature] [L] [0=unordered initial state, 1= ordered] 

Python1.py plots the quantities for L=2, Python2.py for L=20, unordered and ordered initial state. Python3.py plots the histograms. Python4.py plots the quantities as 
function of temperature. 
