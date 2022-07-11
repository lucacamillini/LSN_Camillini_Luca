//
//  esercizio_2.hpp
//  
//
//  Created by luca camillini on 17/03/22.
//

#ifndef esercizio_2_hpp
#define esercizio_2_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "random.h"
#include <fstream>
#include <vector>

//  Parametri per campionamento
unsigned int Camp = 100;

double I;
double I_unif=0.;
double I_2_unif=0.;
double I_quad=0.;
double I_2_quad=0.;
double ave=0.;
double ave_2=0.;
double r[3] = {0,0,0};

Random rnd;

void Eval_unif();
void Eval_quad();
void RW(int k);
void RW_continuum(int k);

#endif /* esercizio_2_hpp */
