#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <armadillo>
#include "time.h"
#include "lib.h"


using namespace std;
using namespace arma;

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
void gauleg(double x1, double x2, double x[], double w[], int N);

int main(int argc, char* argv[])
{


    int N = atoi(argv[1]);
    double *x = new double [N];
    double *w = new double [N];
    double a = -3;
    double b = 3;
    double pi = 3.1415926535;
    double sixteen = 16*16;
    double nevner = 1/sixteen;
    double answer = 5*pi*pi*nevner;
    gauleg(a,b,x,w, N);
    double int_gauss = 0.;
    //   six-double loops
    clock_t start, finish;
    start = clock();

    for (int i=0;i<N;i++){
       for (int j = 0;j<N;j++){
       for (int k = 0;k<N;k++){
       for (int l = 0;l<N;l++){
       for (int m = 0;m<N;m++){
       for (int n = 0;n<N;n++){

            int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);


                }}}}}
        }
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
cout <<  "Running program with N value = " << N << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(10) << "Time used = " << timeused << " seconds." << endl;
cout << setprecision(3) << "The analytical solution is " << int_gauss << endl;
cout << setprecision(3) << "The closed form answer is " <<  answer << endl;
}



double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
    //  this function defines the function to integrate
    double alpha = 2.;
    // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    if (deno < 1e-4){
        return 0;
    }
    else{
        return exp(exp1+exp2)/deno;
    }


} // end of function to evaluate

void gauleg(double x1, double x2, double x[], double w[], int N)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (N + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + N - 1;
   w_low  = w;
   w_high = w + N - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(N + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= N; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = N * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()
