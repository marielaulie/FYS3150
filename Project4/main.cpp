#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <iostream>
#include <armadillo>
#include <fstream>
#include <mpi.h>
using namespace arma;
using namespace std;

// output file
ofstream ofile;

arma::mat make_spinn_matrix(int L);
void initialize_nm_E_M(int L, mat &n_matrix, double &E, double &M);
void Metropolis(int L, int n, mat &n_matrix, double betaT, double &E, double &M,int &acceptance,  vec &ExpectationValue);
void mean(vec &ExpectationValue, int n, int L, double T);

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

int main(int argc, char* argv[])
{
    string filename;
    double T = atof(argv[3]);
    double tempstep = 0.5;
    double initialtemp = 1.0;
    double finaltemp = 3.0;
    int L = atoi(argv[1]); //number of spins
    int n = atoi(argv[2]); //Monte Carlo cycles
    filename=argv[4];
    double M = 0;
    double E = 0;
    double J = 1;
    double betaT = 1.0/T;
    double pi = 3.1415926535;

    string fileout = filename;
    ofile.open(fileout);

    //Defines zero vectors
    vec ExpectationValue = zeros<mat>(5);
    arma::mat n_matrix = make_spinn_matrix(L);

    int acc = 0; //Checking how many flips are accepted

    int numprocs, my_rank;    
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

        initialize_nm_E_M(L, n_matrix, E, M);
        Metropolis(L, n, n_matrix, betaT, E, M,acc, ExpectationValue);
        mean(ExpectationValue, n, L, T);

        std::cout << "Acceptance: " << acc / (double)n /L/L << std::endl;

    MPI_Finalize ();





    //ofile.close();
    return 0;

}


arma::mat make_spinn_matrix(int L){
    //Matrix with the same dimensions as the number of spins
    arma::mat n_matrix = arma::mat(L, L);
    return n_matrix;


}


void initialize_nm_E_M(int L, mat &n_matrix, double &E, double &M){
    //Initializing the spinn matrix with spins in the same direction.
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> distribution(0,1);


    for(int i=0; i < L; i++) {
    for (int j=0; j < L; j++) {
        n_matrix(i,j) = 1.0;
        M +=  (double) n_matrix(i,j);
        //Initializing the spinn matrix with random configuration spins
       // if (n_matrix(i,j) == 0){
         //   n_matrix(i,j) += -1;
       // }

    }
    }
    //initializing the energy and magnetization for the system
    for(int i=0; i < L; i++) {
    for (int j=0; j < L; j++) {

        //multiplies each element in the matrix with the next - using boundary conditions
        E -= (double) n_matrix(i,j)*
                (n_matrix(periodic(i,L,-1),j) +
                 n_matrix(i, periodic(j,L,-1)));



    }
    }

}


void Metropolis(int L, int n, mat &n_matrix, double betaT, double &E, double &M, int &acceptance, vec &ExpectationValue)
{
    //Using the metropolis algorithm to flip spins if the energy is favorable enough

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    arma::vec E_vector = zeros<mat>(17);
    double equilibrium = 5000;


    //Defining the energy vector
    for( int de =-8; de <= 8; de+=4) E_vector(de+8) = exp(-de*betaT) ;
    //for( int de =-8; de <= 8; de++) std::cout << E_vector(de+8) << " " << de << std::endl;


    //Monte Carlo cycles for the Metropolis algorithm
    for (int c = 1; c <= n; c++){
        // loop over all spins
        for(int i =0; i < L; i++) {
        for (int j= 0; j < L; j++){

         // Find random position
         int x = (int) (distribution(gen)*(double)L);
         int y = (int) (distribution(gen)*(double)L);


        //Making an empty vector for the energydiff

        int deltaE =  2*n_matrix(x,y)*
        (n_matrix(x, periodic(y,L,-1))+
         n_matrix(periodic(x,L,-1),y) +
         n_matrix(x, periodic(y,L,1)) +
         n_matrix(periodic(x,L,1),y));
        // Here we perform the Metropolis test



         if ((distribution(gen)) <= E_vector(deltaE+8) ) {
             n_matrix(x,y) *= -1.0;  // flip one spin and accept new spin config
            // update energy and magnetization
             M += (double) 2*n_matrix(x,y);
             E += (double) deltaE;
             acceptance += 1;

      }
   }
}

    //Summing up the values and writing to file

        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << c;
        ofile << setw(15) << setprecision(8) << E/L/L << endl;
        ExpectationValue(0)+=E;
        ExpectationValue(1)+=E*E;
        ExpectationValue(2)+=M;
        ExpectationValue(3)+=M*M;
        ExpectationValue(4)+=fabs(M);

}
    double norm = 1.0/(n);
    ExpectationValue*=norm;

} // end of Metropolis sampling over spins

void mean(vec &ExpectationValue, int n, int L, double T){
    double E_norm = ExpectationValue(0);
    double E2_norm = ExpectationValue(1);
    double M_norm = ExpectationValue(2);
    double M2_norm = ExpectationValue(3);
    double M_abs_norm = ExpectationValue(4);
    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double E_var = (E2_norm - E_norm*E_norm)/L/L;
    double M_var = (M2_norm - M_abs_norm*M_abs_norm)/L/L;
    double kb = 1.;
    double kbT = 1./(kb*T);
    double Cv = E_var*kbT;
    double chi = M_var*kbT;



    cout << "Monte Carlo cycles = " << n << endl;
    cout << "Spinmatrix dimension = " << L << endl;
    cout << "Mean energy = " << E_norm/L/L << endl;
    cout << "Mean magnetization = " << M_norm << endl;
    cout << "Mean abs magn = " << M_abs_norm << endl;
    cout << "Energy variance = " << E_var << endl;
    cout << "Magnetization variance = " << M_var << endl;
    cout << "The heat capasity of the system = " << Cv << endl;
    cout << "The susceptibility of the system = " << chi << endl;




}
