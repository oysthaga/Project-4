# include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <cassert>

int PBC(int i, int L)
{
    return (i+L)%L;
}


int main(int argc, char* argv[])
{

    
    std::string in = argv[1];
    int NumCycles = pow(10,atoi(argv[1])); // Number of Monte Carlo cycles
    double T = atof(argv[2]); // Temperature
    double beta = 1./T;
    int L = atoi(argv[3]); // Length of grid 
    int N = L*L;
    int in4 = atoi(argv[4]); // Ordered or unordered initial condition
    arma::imat s;
    if (in4==0) // Unordered initial state
    {
        s = arma::randi(L,L, arma::distr_param(0,1)); // Entries random 0 or 1
        s.replace(0,-1); // Replace 0 with -1
    }
    else if (in4==1) // Ordered initial state
    {
        s.ones(L,L); // All ones 
    }
    else
    {
        throw std::invalid_argument( "4th input must be 0 or 1." );
    }

    unsigned int seed = 1432890;
    std::mt19937 generator; // Random value genenerator
    generator.seed(seed);
    std::uniform_int_distribution<int> uniform1(0,L-1);
    std::uniform_real_distribution<double> uniform2(0.0,1.0);

    // Array for p-ratio so we don't have to call exp in the for-loop. 
    arma::vec exparray = arma::zeros(16+1);
    exparray(0) = exp(-beta*(-8));
    exparray(4) = exp(-beta*(-4));
    exparray(8) = 1;
    exparray(12) = exp(-beta*4);
    exparray(16) = exp(-beta*8);

    arma::vec epsilon = arma::vec(NumCycles); 
    arma::vec m = arma::vec(NumCycles); 
    arma::vec epsilon_squared = arma::vec(NumCycles); 
    arma::vec m_squared = arma::vec(NumCycles); 
    
    double E0 = 0;
    arma::vec E = arma::vec(NumCycles+1);
    for (int i=0; i < L; i++)
    {
        for (int j=0; j < L; j++)
        {
            E0 -=s(i,j)*( s(i,PBC(j+1, L))
                        + s(PBC(i+1, L),j));
        }
    }
    E(0) = E0;


    
    for (int j = 0; j < NumCycles; j++)
    {
        E(j+1) = E(j);
        for (int i = 0; i <= N; i++)
        {
            int n = uniform1(generator);
            int m = uniform1(generator);

           
    
            arma::imat s_prime = s;
            s_prime(n, m) = -s(n, m);

           
    
            double Ei        = -s(n,m)*(   s(n,PBC(m-1, L)) 
                                        + s(n,PBC(m+1, L))
                                        + s(PBC(n-1, L),m)
                                        + s(PBC(n+1, L),m));

            double Ei_prime  = -s_prime(n,m)*(  s(n,PBC(m-1, L)) 
                                            + s(n,PBC(m+1, L))
                                            + s(PBC(n-1, L),m)
                                            + s(PBC(n+1, L),m));
            
            double DeltaE = Ei_prime -Ei;
            double p_ratio = exparray(DeltaE+8);
            
            
            double A = std::min(1., p_ratio);
            double r = uniform2(generator);
            if (r < A)
            {
                s = s_prime;
                E(j+1) += DeltaE;
            }
        }
        epsilon = E/N;
        m(j) = arma::accu(s)/(double)N;
        
        
    }

    double mean_epsilon = arma::mean(epsilon);
    double mean_epsilon_squared = arma::mean(epsilon%epsilon);
    double mean_abs_m = arma::mean(arma::abs(m));
    double mean_m_squared = arma::mean(m%m);
    
    double C_V = (N/(pow(T,2)))
                *(    mean_epsilon_squared
                    - mean_epsilon*mean_epsilon); // [k_B]
    double chi = (N/T)
                *(    mean_m_squared
                    - mean_abs_m*mean_abs_m); // [1/J]


    std::cout << "mean_epsilon";
    std::cout << '\n';
    std::cout << mean_epsilon;
    std::cout << '\n';
    std::cout << "mean_epsilon_squared";
    std::cout << '\n';
    std::cout << mean_epsilon_squared;
    std::cout << '\n';
    std::cout << "mean_abs_m";
    std::cout << '\n';
    std::cout << mean_abs_m;
    std::cout << '\n';
    std::cout << "mean_m_squared";
    std::cout << '\n';
    std::cout << mean_m_squared;
    std::cout << '\n';
    std::cout << "C_V";
    std::cout << '\n';
    std::cout << C_V;
    std::cout << '\n';
    std::cout << "chi";
    std::cout << '\n';
    std::cout << chi;
    std::cout << '\n';
    
    std::string name = "FileMat";
    name += "in";
    name += in;
    name += "T";
    name += argv[2];
    name += "L";
    name += argv[3];
    name += "ord";
    name += argv[4];
    name += ".bin";
    arma::vec FileMat = arma::vec{mean_epsilon, mean_epsilon_squared, mean_abs_m, mean_m_squared, C_V, chi};
    FileMat.save(name); 


    std::string name1 = "epsilon";
    name1 += "in";
    name1 += in;
    name1 += "T";
    name1 += argv[2];
    name1 += "L";
    name1 += argv[3];
    name1 += "ord";
    name1 += argv[4];
    name1 += ".bin";
    epsilon.save(name1); 


    std::string name2 = "m";
    name2 += "in";
    name2 += in;
    name2 += "T";
    name2 += argv[2];
    name2 += "L";
    name2 += argv[3];
    name2 += "ord";
    name2 += argv[4];
    name2 += ".bin";
    m.save(name2); 
    

}