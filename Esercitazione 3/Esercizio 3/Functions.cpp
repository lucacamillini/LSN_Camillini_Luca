//
//  Functions.cpp
//  
//
//  Created by luca camillini on 31/03/22.
//

#include <stdio.h>
#include "Esercizio3.hpp"

double call(double S, double T, double K, double sigma, double r){
    double d_1, d_2, call;
    d_1 = 1/(sigma*sqrt(T)) * (log(S/K) + (r + pow(sigma,2)/2 * T));
    d_2 = d_1 - sigma*sqrt(T);
    
    call = S*(0.5 * (1 + erf(d_1/sqrt(2)))) - K*exp(-r*T)*(0.5 * (1 + erf(d_2/sqrt(2))));
    return call;
}

double put(double S, double T, double K, double sigma, double r){
    double d_1, d_2, pull;
    d_1 = 1/(sigma*sqrt(T)) * (log(S/K) + (r + pow(sigma,2)/2 * T));
    d_2 = d_1 - sigma*sqrt(T);
    
    pull = S*(0.5 * (1 + erf(d_1/sqrt(2))) - 1) - K*exp(-r*T)*(0.5*(1 + erf(d_2/sqrt(2))) - 1);
    
    return pull;
}


