# include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <cassert>

int PBC(int i, int L) // Periodic boundary conditions
{
    return (i+L)%L;
}


int main(int argc, char* argv[])
{

    
    std::string in = argv[1]; 
    int NumCycles = pow(10,atoi(argv[1])); // Number of Monte Carlo cycles
    int L = atoi(argv[2]); // Length of grid
    int N = L*L;
    int in4 = atoi(argv[3]); // Ordered or unordered initial condition
    int K = atoi(argv[4]); // Number of temperature steps
    double Tmin = atof(argv[5]); // Minimum temperature value
    double Tmax = atof(argv[6]); // Maximum temperature value 

    unsigned int seed = 1432890;
    std::mt19937 generator; // Random value genenerator
    generator.seed(seed);
    std::uniform_int_distribution<int> uniform1(0,L-1);
    std::uniform_real_distribution<double> uniform2(0.0,1.0);
    

    arma::vec T = arma::linspace(Tmin, Tmax, K+1);
    arma::vec beta = 1./T;

    arma::vec mean_epsilon = arma::vec(K+1);
    arma::vec mean_epsilon_squared = arma::vec(K+1);
    arma::vec mean_abs_m = arma::vec(K+1);
    arma::vec mean_m_squared = arma::vec(K+1);
    
    arma::vec C_V = arma::vec(K+1);
    arma::vec chi = arma::vec(K+1);

    // Parallellization
    #ifdef _OPENMP
    {
        #pragma omp parallel for
        for (int k = 0; k <= K; k++)
        {
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

            // Initial energy
            double E0 = 0;
            for (int i=0; i < L; i++)
            {
                for (int j=0; j < L; j++)
                {
                    E0 -=s(i,j)*( s(i,PBC(j+1, L))
                                + s(PBC(i+1, L),j));
                }
            }

            arma::vec epsilon = arma::vec(NumCycles); 
            arma::vec m = arma::vec(NumCycles); 

            arma::vec E = arma::vec(NumCycles+1);
            E(0) = E0;

            // Array for p-ratio so we don't have to call exp in j-loop. 
            arma::vec exparray = arma::zeros(16+1);
            exparray(0) = exp(-beta(k)*(-8));
            exparray(4) = exp(-beta(k)*(-4));
            exparray(8) = 1;
            exparray(12) = exp(-beta(k)*4);
            exparray(16) = exp(-beta(k)*8);
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

            mean_epsilon(k) = arma::mean(epsilon);
            mean_epsilon_squared(k) = arma::mean(epsilon%epsilon);
            mean_abs_m(k) = arma::mean(arma::abs(m));
            mean_m_squared(k) = arma::mean(m%m);
            
            C_V(k) = (N/(pow(T(k),2)))
                        *(    mean_epsilon_squared(k)
                            - mean_epsilon(k)*mean_epsilon(k)); // [k_B]
            chi(k) = (N/T(k))
                        *(    mean_m_squared(k)
                            - mean_abs_m(k)*mean_abs_m(k)); // [1/J]

           
            


        }
    }
    #else
    {
        // Identical, this part is compiled if we do not parallelize.
        for (int k = 0; k <= K; k++)
        {
            arma::vec E = arma::vec(NumCycles+1);
            E(0) = E0;
            arma::vec exparray = arma::zeros(16+1);
            exparray(0) = exp(-beta(k)*(-8));
            exparray(4) = exp(-beta(k)*(-4));
            exparray(8) = 1;
            exparray(12) = exp(-beta(k)*4);
            exparray(16) = exp(-beta(k)*8);
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

            mean_epsilon(k) = arma::mean(epsilon);
            mean_epsilon_squared(k) = arma::mean(epsilon%epsilon);
            mean_abs_m(k) = arma::mean(arma::abs(m));
            mean_m_squared(k) = arma::mean(m%m);
            
            C_V(k) = (N/(pow(T(k),2)))
                        *(    mean_epsilon_squared(k)
                            - mean_epsilon(k)*mean_epsilon(k)); // [k_B]
            chi(k) = (N/T(k))
                        *(    mean_m_squared(k)
                            - mean_abs_m(k)*mean_abs_m(k)); // [1/J]

        }
    }
    #endif

    std::string name1 = "epsilon";
    name1 += "in";
    name1 += in;
    name1 += "L";
    name1 += argv[2];
    name1 += "ord";
    name1 += argv[3];
    name1 += "Tmin";
    name1 += argv[5];
    name1 += "Tmax";
    name1 += argv[6];
    name1 += ".bin";
    mean_epsilon.save(name1); 


    std::string name2 = "epsilon_squared";
    name2 += "in";
    name2 += in;
    name2 += "L";
    name2 += argv[2];
    name2 += "ord";
    name2 += argv[3];
    name2 += "Tmin";
    name2 += argv[5];
    name2 += "Tmax";
    name2 += argv[6];
    name2 += ".bin";
    mean_epsilon_squared.save(name2); 

    std::string name3 = "m";
    name3 += "in";
    name3 += in;
    name3 += "L";
    name3 += argv[2];
    name3 += "ord";
    name3 += argv[3];
    name3 += "Tmin";
    name3 += argv[5];
    name3 += "Tmax";
    name3 += argv[6];
    name3 += ".bin";
    mean_abs_m.save(name3); 


    std::string name4 = "m_squared";
    name4 += "in";
    name4 += in;
    name4 += "L";
    name4 += argv[2];
    name4 += "ord";
    name4 += argv[3];
    name4 += "Tmin";
    name4 += argv[5];
    name4 += "Tmax";
    name4 += argv[6];
    name4 += ".bin";
    mean_m_squared.save(name4); 

    std::string name5 = "C_V";
    name5 += "in";
    name5 += in;
    name5 += "L";
    name5 += argv[2];
    name5 += "ord";
    name5 += argv[3];
    name5 += "Tmin";
    name5 += argv[5];
    name5 += "Tmax";
    name5 += argv[6];
    name5 += ".bin";
    C_V.save(name5); 

    std::string name6 = "chi";
    name6 += "in";
    name6 += in;
    name6 += "L";
    name6 += argv[2];
    name6 += "ord";
    name6 += argv[3];
    name6 += "Tmin";
    name6 += argv[5];
    name6 += "Tmax";
    name6 += argv[6];
    name6 += ".bin";
    chi.save(name6); 

    std::string name7 = "T";
    name7 += "k";
    name7 += argv[4];
    name7 += "Tmin";
    name7 += argv[5];
    name7 += "Tmax";
    name7 += argv[6];
    name7 += ".bin";
    T.save(name7); 


    
}