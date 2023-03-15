/* Generating the 1/f^a family of noise*/

#include <random>
#include <cmath>
#include "noise.h" //comment out for python interfacing
#include <vector>
#include <iostream>

extern "C"  //works with ctypes

/*

addpoles() adds poles contingent on the function check().  
For now, adds 3 poles per decade above 10^6 steps.

*/

signed long addpoles(signed long p, double d){

p = ceil(p + d);

return (signed long)p;
}

/*

filters() is a 2d vector which stores the filter coefficients.
Expands by n poles, contingent on check()==TRUE.

*/

std::vector<std::vector<double>> filters(signed long poles, double gam, double ca){


std::vector<std::vector<double>> vec (2, std::vector<double>(poles));

for(signed long i = 0; i<2; i++){

if(i == 0){

for(signed long j = 0; j<poles; j++){

vec[i][j] = exp(-2.0*M_PI*pow(10.0, ((0.5*(j - poles)) - (gam/4.0) - ca)));

}
}

if(i == 1){

for(signed long j = 0; j<poles; j++){

vec[i][j] = exp(-2.0*M_PI*pow(10.0, ((0.5*(j - poles)) - ca)));

}

}

}


return vec;

}

/*

check() is a fairly fast method of checking if the steps of a signal
has reached a new decade. Much faster than log10(steps) mod 1 == 0.
Will check up to one quintillion -- should be sufficient for most applications.

*/

bool check(signed long s){

bool bt = true;
bool bf = false;

if(
s == 1000000 ||
s == 10000000 ||
s == 100000000 ||
s == 1000000000 ||
s == 10000000000 ||
s == 100000000000 ||
s == 1000000000000 ||
s == 10000000000000 ||
s == 100000000000000 ||
s == 1000000000000000 ||
s == 10000000000000000 ||
s == 100000000000000000 ||
s == 1000000000000000000    
){

return bt;

} else {

return bf;

}
}

/* 

Main function begins here. 

steps = desired length of a signal

gamma = alpha [0, 2]  (often interchanged)

fhigh = sampling frequency 

flow = lowest frequency 

Gain = scalar of noise

h = number of poles per decade

*/

std::vector<double> noise(signed long steps, double gamma, double fhigh, double flow, double Gain, double h){	
std::random_device rd;                                    //random_device is a seed generator
std::default_random_engine generator(rd());               //source of randomness
std::normal_distribution<double> distribution(0.0, 1.0);  //0 mean and unity sd, rnorm

//DEFAULT VALUES 

double DC = 0.0;
double corneradjust = 0.0;
double m = log10(fhigh/flow); 		  //number of decades
signed long n = ceil(h*m);            //total number of poles
signed long f = 3;              	  //filter n times. "Primes the pump" - E.V.
signed long i, j, k;

std::vector<double> x(n+1), y(n), xold(n), yold(n); //a, b used to be here
std::vector<std::vector<double>> fil;

fil = filters(n, gamma, corneradjust);

x.reserve(n+1);
y.reserve(n);
xold.reserve(n);
yold.reserve(n);

std::vector<double> NoiseOut;
NoiseOut.reserve(steps);

for (i=0; i<n; i++){
xold[i] = 0.0;
yold[i] = 0.0;
}

for (j=0;j<f*n;j++){

x[0] = distribution(generator);

for (i=0;i<n;i++){
y[i] = x[i] + fil[0][i]*yold[i] - fil[1][i]*xold[i];
x[i+1] = y[i];
xold[i] = x[i];
yold[i] = y[i];
}
}

for (k=0;k<steps; k++){ 

if(check(k) == true){

n = addpoles(n, h);
fil = filters(n, gamma, corneradjust);

x.insert(x.begin(), h, 0);
xold.insert(xold.begin(), h, 0);
y.insert(y.begin(), h, 0);
yold.insert(yold.begin(), h, 0);

y.reserve(n);
yold.reserve(n);
x.reserve(n+1);
xold.reserve(n);

}

x[0] = distribution(generator); 

for (i=0;i<n;i++){
y[i] = x[i] + fil[0][i]*yold[i] - fil[1][i]*xold[i];
x[i+1] = y[i];
xold[i] = x[i];
yold[i] = y[i];
}

NoiseOut[k] = DC + Gain*y[n-1];

}

return NoiseOut;

}
