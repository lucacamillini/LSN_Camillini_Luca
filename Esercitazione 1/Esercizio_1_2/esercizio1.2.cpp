//
//  esercizio1_2.cpp
//  
//
//  Created by luca camillini on 24/03/22.
//

#include "esercizio1.2.hpp"

int main(){
    
    //parametri
    int realizations = 10000;
    int samplings[4] = {1,2,10,100};

    //inizializzo il generatore
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
    
//    apro il file di scrittura
    ofstream out;
    out.open("1_sampl.csv");
    
    for (int i=0; i<realizations; i++) {
        int N = 1;
        double S_unif=0, S_exp=0, S_lorentz=0;
        
        for (int j=0; j<N; j++) {
            S_unif += rnd.Rannyu(1,6);
        }
        for (int j=0; j<N; j++) {
            S_exp += rnd.Exp(1);
        }
        for (int j=0; j<N; j++) {
            S_lorentz += rnd.Lorentz(0.,1.);
        }
        out << S_unif/N << ", " << S_exp/N << ", " << S_lorentz/N << endl;
    }
    out.close();
    
//    apro il file di scrittura
    out.open("2_sampl.csv");
    
    for (int i=0; i<realizations; i++) {
        int N = 2;
        double S_unif=0, S_exp=0, S_lorentz=0;
        
        for (int j=0; j<N; j++) {
            S_unif += rnd.Rannyu(1,6);
            S_exp += rnd.Exp(1);
            S_lorentz += rnd.Lorentz(0,1);
        }
        out << S_unif/N << ", " << S_exp/N << ", " << S_lorentz/N << endl;
    }
    out.close();

//    apro il file di scrittura
    out.open("10_sampl.csv");
    
    for (int i=0; i<realizations; i++) {
        int N = 10;
        double S_unif=0, S_exp=0, S_lorentz=0;
        
        for (int j=0; j<N; j++) {
            S_unif += rnd.Rannyu(1,6);
            S_exp += rnd.Exp(1);
            S_lorentz += rnd.Lorentz(0,1);
        }
        out << S_unif/N << ", " << S_exp/N << ", " << S_lorentz/N << endl;
    }
    out.close();

//    apro il file di scrittura
    out.open("100_sampl.csv");
    
    for (int i=0; i<realizations; i++) {
        int N = 100;
        double S_unif=0, S_exp=0, S_lorentz=0;
        
        for (int j=0; j<N; j++) {
            S_unif += rnd.Rannyu(1,6);
            S_exp += rnd.Exp(1);
            S_lorentz += rnd.Lorentz(0,1);
        }
        out << S_unif/N << ", " << S_exp/N << ", " << S_lorentz/N << endl;
    }
    out.close();
    
    out.open("distribution_example.csv");
    for (int i=0; i<10000; i++) {
        double S_unif=0, S_exp=0, S_lorentz=0;
        S_unif = rnd.Rannyu(1,6);
        S_exp = rnd.Exp(1);
        S_lorentz = rnd.Lorentz(0,1);
        out << S_unif << ", " << S_exp << ", " << S_lorentz << endl;
    }
    out.close();
    
    return 0;
}
