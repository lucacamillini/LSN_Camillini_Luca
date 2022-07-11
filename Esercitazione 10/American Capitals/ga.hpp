//
//  ga.hpp
//  
//
//  Created by luca camillini on 23/05/22.
//

#ifndef ga_hpp
#define ga_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include "random.h"

using namespace std;

class ga{
    
private:
    int m_N, m_type, m_pop; //Number of city, circle or square lattice, popolation size.
    double m_p_cross, m_p_mut;
    Random m_rnd; //random generator
    vector<vector<double> > m_world; //cities configuration
    vector<vector<int> > m_population; //population

public:
    //Constructor
    ga(int, int, int, Random, double, double); //Prepares the cities in the region prescribed. The first term is the number of cities, with the second int I can choose between the circle or the square. 0 = circle, 1 = square. The third int is the size of the population. Last is the random number generator. The construcor also creates the initial population. Doubles are the crossover and mutations probabilities.
    //Destructure
    ~ga(){;};

    void mutation(vector<int>&);
    int selection();
    void Swap(int);
    double cost_function(vector<int>);
    
    void check(vector<int>&);
    
    void First_Gen();
    void Sort();
    void Next_Gen();
    
    vector< vector<int> > Get_m_pop();
    vector< vector<double> > Get_m_world();
    vector<int> Get_Chrom(int);
    
    void Set_Chrom(vector<int>, int);
    void Set_World(string a);
    
    void print(string , string );
    void print_world(string);
    void print_half(string);
    
};

#endif /* ga_hpp */
/*Che mutazioni si possono fare? Swap, riscrittura del gene a specchio, invertire il prima col dopo mantenendo la posizione 1 al suo posto, individuare due estremi e fare lo swap degli estremi mentre nelle parti in mezzo invertire l'ordine*/
