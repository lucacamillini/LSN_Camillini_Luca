//
//  main.cpp
//  
//
//  Created by luca camillini on 24/05/22.
//

#include "main.hpp"

using namespace std;

int main(){
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    cout << "=============Salesman Problem and GENETIC ALGORITHM=============" << endl;
    
    ga prova(34,1,200, rnd, 0.5, 0.1);
    
    prova.Sort();
    prova.print();

    for(int i=0; i<1000; i++){
        if((i+1)%10==0) cout << "Generazione " << i+1 << endl;
        prova.Next_Gen();
        prova.print();
        prova.print_half();
    }
    return 0;
}
