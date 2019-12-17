#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <iostream>
#include <armadillo>
#include <fstream>
using namespace arma;
using namespace std;

arma::mat make_BE(int n);
void initialize_BE(int n, mat &BE, vec &b, vec &a, vec &u, double alpha);
void gaus_eliminate(int n, vec diag, vec offdiag,  double alpha);
void forward_step(int n, double alpha, vec &v, vec &vprev);
//void forward_step(int n, vec u, vec uprev, double alpha);

int main(int argc, char* argv[]){
    cout << "Hello World!" << endl;

    int T = 100;
    double n = 9;
    double alpha = 1.0;//tstep/(xstep*xstep);
    double x = 0.0;
    vec b = zeros<mat>(n);
    vec u = zeros<mat>(n+2);
    vec a = zeros<mat>(n-1);
    vec v= zeros<vec>(n);
    vec vprev= zeros<vec>(n);
    arma::mat BE = make_BE(n);


    //gaus_eliminate(n, diag, offdiag, alpha);


    //backward euler


        initialize_BE(n, BE, b, a, u, alpha);
        forward_step(n, alpha, v, vprev);



}

void forward_step(int n, double alpha, vec &v, vec &vprev){
    //for (int j=1; j < n; j++){
    for (int i = 0; i < n; i++){
        vprev(i)=3.0;
        v(i) = alpha*vprev(i-1)+(1-2*alpha)*vprev(i)+alpha*vprev(i+1);
    }
//}
    cout << v << endl;
}

arma::mat make_BE(int n){
    //Zero matrix
    arma::mat BE = arma::mat(n, n);
    return BE;

}



void initialize_BE(int n, mat &BE, vec &b, vec &a, vec &u, double alpha){
    double diag = 1+2*alpha;
    double offdiag = -alpha;

    //initializing the vector with the diagonal elements
    for (int i=0; i < n; i++){
        b(i) = diag;
    }
    //initializing the vector with the super and sub-diagonal elements.
    for (int i=0; i < n-1; i++){
        a(i) = offdiag;
    }
    //testing if it works
    cout << diag << endl;
    cout << offdiag << endl;


    for (int i = 1; i < n; i++){
        a(i-1) = a(i-1)/b(i-1);
        u(i) = u(i)/b(i-1);
        b(i-1) = 1.0;

        u(i+1) += u(i)*alpha;
        b(i) += a(i-1)*alpha;
        u(n)= u(n)/b(n-1);
        b(n-1)=1.0;

    }
    for (int i = n; i > 1; i--){
        u(i-1) -=u(i)*a(i-2);


    }
    cout << b << endl;
    cout << a << endl;
    cout << u << endl;



}

