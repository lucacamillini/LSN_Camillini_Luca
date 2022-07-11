//
//  Esercizio3.hpp
//  
//
//  Created by luca camillini on 31/03/22.
//

#ifndef Esercizio3_hpp
#define Esercizio3_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include <vector>

using namespace std;

// Blocking parameters
int M = 1e5;
int N = 200, L=M/N;
double C_sum_prog=0, C_sum_prog_2=0, P_sum_prog=0, P_sum_prog_2=0, C_blk_ave=0, P_blk_ave=0;

//Random generator
Random rnd;

//Parameters
double S_0 = 100;
double K = 100;
double T = 1;
double r = 0.1;
double sigma = 0.25;

//Files
ofstream out_1, out_2, out_3, out_4;

void Inizialization();
void Reset(int);
double GBM(int n_steps, Random rnd, double mu, double sigma);

#endif /* Esercizio3_hpp */
