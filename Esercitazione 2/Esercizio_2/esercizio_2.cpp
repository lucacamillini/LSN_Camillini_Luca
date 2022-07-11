//
//  esercizio_2.cpp
//  
//
//  Created by luca camillini on 17/03/22.
//

#include "esercizio_2.hpp"

using namespace std;

int main(){
    
//  Parametri per data blocking
    int M = 100000;
    int N = 150;
    int L = M/N;
    
//    random generator
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
    
//  output files
    ofstream out_1, out_3;
    out_1.open("integrale_media.csv");
    out_3.open("integrale_quad.csv");

//  02.1
//  uniform distribution sampling
    for(int j = 1; j<=N; j++){
        I=0;
        for(int i=0; i<L; i++){
            Eval_unif();
        }
        I_unif += I/(L*Camp);
        I_2_unif += pow(I/(L*Camp),2);
        if(j==1){
            out_1 << j << ", " << I_unif/double(j) << ", " << 0 << endl;
        }
        else{
            out_1 << j << ", " << I_unif/double(j) << ", " << sqrt(abs(I_2_unif/double(j) - pow(I_unif/double(j),2))/double(j-1)) << endl;
        }
    }
    
//  importance sampling
//  Here it's used p(x) = 3/2 * (1 - x^2) and g(x) = (M_PI/3) * cos((M_PI/2.) * x) * 1./(1 - x^2);

//  Inizializzo le somme progressive
    I_unif = 0;
    I_2_unif = 0;
    
    for(int j = 1; j<=N; j++){
        I=0;
        for(int i=0; i<L; i++){
            Eval_quad();
        }
        I_unif += I/(L*Camp);
        I_2_unif += pow(I/(L*Camp),2);
        if(j==1){
            out_3 << j << ", " << I_unif/double(j) << ", " << 0 << endl;
        }
        else{
            out_3 << j << ", " << I_unif/double(j) << ", " << sqrt(abs(I_2_unif/double(j) - pow(I_unif/double(j),2))/double(j-1)) << endl;
        }
    }
    
    out_1.close();
    out_3.close();
    
//  02.2.1
    
//  output files
    ofstream out_2;
    out_2.open("cubic_lattice.csv");
    cout << "RW" << endl;
//  Lunghezza RW e parametri per blocking
    int n = 100;
    M = 10000;
    N = 100;
    L = M/N;

    for (int k = 0; k<n; k++) {
        ave = 0;
        ave_2 = 0;
        for (int i=1; i<=N; i++) {
            double rw = 0;
            for(int h=0; h<L; h++){
                double r_2 = 0;
                r[0] = 0;
                r[1] = 0;
                r[2] = 0;
                RW(k);
                rw += sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));
            }
            ave += rw/L;
            ave_2 += pow(rw/L,2);
        }
        out_2 << k << ", " << ave/N << ", " << sqrt(abs(ave_2/N - pow(ave/N,2))/(N-1)) << endl;
    }
    out_2.close();
    
//  02.2.2
     M = 10000;
     N = 100;
     L = M/N;
 
    ofstream out_4;
    out_4.open("continous_lattice.csv");
    
    for (int i=1; i<=n ; i++) {
        ave = 0;
        ave_2 = 0;
//      incomincio il ciclo per mediare
        for (int j=0; j<N; j++) {
            double rw = 0;
            for(int h=0; h<L; h++){
                r[0] = 0;
                r[1] = 0;
                r[2] = 0;
                RW_continuum(i);
                rw += sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));
            }
            ave += rw/L;
            ave_2 += pow(rw/L,2);
        }
        out_4 << i << ", " << ave/N << ", " << sqrt(abs(ave_2/N - pow(ave/N,2))/(N-1)) << endl;
    }
    
    out_4.close();

    return 0;
}

void Eval_unif(){
    for(int i=0; i<Camp; i++){
        I += (M_PI/2.) * cos(M_PI*rnd.Rannyu()/2.);
    }
}

void Eval_quad(){
    for(int i=0; i<Camp; i++){
        double x=rnd.Quad();
        I += (M_PI/3.) * cos((M_PI/2.)*x)*1./(1-pow(x,2));
    }
}

void RW(int k){
    for (int j=0; j<k; j++) {
        int index = int(rnd.Rannyu(0,3));
        float sign;
        sign = rnd.Rannyu(-1,1);
        sign = int(sign/abs(sign));
        r[index] += sign;
    }
}

void RW_continuum(int i){
    for (int k=0; k<i; k++) {
        double theta, phi;
        theta = rnd.Rannyu(0,M_PI);
        phi = rnd.Rannyu(0,2.*M_PI);
        r[0] += sin(theta)*cos(phi);
        r[1] += sin(theta)*sin(phi);
        r[2] += cos(theta);
    }
}
