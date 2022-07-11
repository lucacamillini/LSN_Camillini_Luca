//
//  esercizio1.cpp
//  
//
//  Created by luca camillini on 24/03/22.
//

#include "esercizio1.3.hpp"

int main(){
    
    //parametri
    int M = 1000000;
    int N = 100;
    int L = M/N;

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
    
    double sum_ave=0, sum_ave_2=0;
    
    ofstream out;
    out.open("pi_value.csv");
    
    for (int i=0; i<N; i++) {
        double l=0.875;
        int counter = 0;
        for (int j=0; j<L; j++) {
            double y_0=rnd.Rannyu();
            double x,y;
            x = rnd.Rannyu(-1,1);
            y = rnd.Rannyu();
            while (x*x + y*y > 1) {
                x = rnd.Rannyu(-1,1);
                y = rnd.Rannyu();
            }
            double sin_theta = y/sqrt(x*x + y*y);
            if (y_0 + 0.5*l*sin_theta >= 1 or y_0 - 0.5*l*sin_theta <= 0) {
                counter += 1;
            }
        }
        sum_ave += 2*l*double(L)/double(counter);
        sum_ave_2 += pow(2*l*L/counter,2);
        if (i==0) {
            out << sum_ave/i << ", " << 0 << endl;
        }
        else{
            out << sum_ave/(i+1) << ", " << sqrt(abs(sum_ave_2/(i+1) - pow(sum_ave/(i+1),2))/i) << endl;
        }
    }
    
    out.close();
    return 0;
}
