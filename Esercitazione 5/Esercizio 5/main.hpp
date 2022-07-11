//
//  main.hpp
//  
//
//  Created by luca camillini on 28/04/22.
//

//#ifndef main_hpp
//#define main_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

using namespace std;

//Parameters
double a = 0.0529;
double step = 3*a;
unsigned int nstep = 100, i = 0, counter = 0;
double r_old[3] = {1,1,1};
double r_new[3];

Random rnd;

//Blocking
unsigned int M = 1000000;
unsigned int N = 100;
unsigned int L = M/N;
double sum_prog=0, sum_prog_2=0, r_0=0, r_mean=0;

//Probabilities
double A;

//Files
ofstream out_1, out_2, out_3, out_4;

//Functions
void Evolution(unsigned int d);
void passo();
double min(double a, double b);
void A_0(unsigned d);
double p_1(double h);
double p_2(double h, double r_z);
void write_conf();
void Inizialization();
void Reset();


//#endif /* main_hpp */
