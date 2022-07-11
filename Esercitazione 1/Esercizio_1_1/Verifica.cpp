//
//  Verifica.cpp
//  
//
//  Created by luca camillini on 16/03/22.
//

#include <stdio.h>
#include "esercizio1.h"

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
    
    vector<double> dati;
    int M=10000;
    
    for(int i=0; i<M; i++){
        dati.push_back(rnd.Lorentz(5,2));
    }
    
    ofstream out;
    out.open("dati_prova.csv");
    for(int i=0; i<dati.size(); i++){
        out << dati[i] << endl;
    }
    
    out.close();
    rnd.SaveSeed();
    
    return 0;
}
