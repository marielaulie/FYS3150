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
#include <math.h>
using namespace arma;
using namespace std;

arma::mat make_u(int n);
arma::mat make_utemp(int n);
void initialize_uu_xy(mat u, mat utemp, int n, vec x, vec y, double h);
void Jabobi_Laplace(mat u, mat utemp, int n);
void exact(int n, vec x, vec y, double h);
void jacobi(int n, mat u, mat utemp, double h, vec x, vec y);

int main(int argc, char* argv[])
{
    int n = atof(argv[1]);
    double h = 1/double((n+1));
    double L = 1.0;
    vec x = zeros<mat>(n+1);
    vec y = zeros<mat>(n+1);
    arma::mat u = make_u(n);
    arma::mat utemp = make_utemp(n);
    initialize_uu_xy(u, utemp, n, x, y, h);
    //Jabobi_Laplace(u, utemp, n);
    exact(n, x, y,h);
    jacobi(n, u, utemp, h, x, y);



    return 0;
}
arma::mat make_u(int n){
    //Matrix with the same dimensions as the number of spins
    arma::mat u = arma::mat(n+1, n+1);

    return u;


}
arma::mat make_utemp(int n){
    //Matrix with the same dimensions as the number of spins
    arma::mat utemp = arma::mat(n+1, n+1);

    return utemp;


}

void initialize_uu_xy(mat u, mat utemp, int n, vec x, vec y, double h){
    //implementing boundary condtitions u(0,y,t)=u(L,y,t)=u(x,0,t)=u(x,L,t)=0
    for (int i = 1; i < n; i++){
        for (int j = 1; j < n; j++){
            utemp(i,j) = 1.0;


        }
    }

    for (int i = 0; i < n+1; i++){
        for (int j = 0; j < n+1; j++){
            u(i,0)=0;
            u(i, n)=0;
            u(0,j) = 0;
            u(n,j) = 0;


        }
    }

    cout << u << endl;
}

void exact(int n, vec x, vec y, double h){
    double exact_solution;
    double pi = 3.141595;
    for (double i = 0; i< n+1; i++){
        x(i) = i*h;
        y(i) = x(i);
        exact_solution =  -double (sin(x(i)*pi))*double(sin(y(i)*pi));
        cout << exact_solution << endl;

    }
    cout << exact_solution << endl;


}


void jacobi(int n, mat u, mat utemp, double h, vec x, vec y){
    double dt = 0.25*h*h;
    double pi = 3.141595;
    double alpha = dt/h*h;
    arma::mat q = arma::mat(n+1, n+1);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            x(i) = i*h;
            y(i) = x(i);
            q(i,j) = -double (sin(x(i)*pi))*double(sin(y(i)*pi));

        }
    }
    for (int t = 0; t < 1430; t++){
        for (int i = 1; i < n; i++){
          for (int j = 1; j < n; j++){
              u(i,j) = double(dt*q(i,j)) + utemp(i,j) + alpha*(utemp(i+1, j) + utemp(i-1,j) + utemp(i,j+1) + utemp(i,j-1)-4*utemp(i,j));

         }
      }

        double sum = 0.0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                sum += (utemp(i,j)-u(i,j))*(utemp(i,j)-u(i,j));
                utemp(i,j) = u(i,j);


}
        }


        if(sqrt (sum) < 0.00001){
            cout <<  t << "=" << u << endl;

        }


}








}




/*
void Jabobi_Laplace(mat u, mat utemp, int n){
    int iterations = 0;
    int max_it = 10000000;
    double converge_tol = 0.00000001;
    double diff;
    while ((iterations <= max_it) && (diff > converge_tol )){
           diff = 0.0;
           utemp = u;
           for (int j = 1; j <= n; j++){
               for (int l = 1; l <=n; l++){
                   u(j,l) = 0.25*(utemp(j+1, l) + utemp(j-1,l) + utemp(j,l+1) + utemp(j,l-1));
                   diff += fabs(utemp(j,l)-u(j,l));


               }
           }
           iterations += 1;
           diff/=pow(n,2.0);
}

cout << iterations << endl;

}*/
