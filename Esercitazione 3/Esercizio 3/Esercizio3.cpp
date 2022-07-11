//
//  Esercizio3.cpp
//  
//
//  Created by luca camillini on 31/03/22.
//

#include "Esercizio3.hpp"

int main(){
//  inizializzo il generatore di numeri casuali
    Inizialization();
    
    cout << "================== Plain vanilla option pricing ==================" << endl;
    
    out_1.open("call.csv");
    out_2.open("put.csv");
    out_3.open("C_simulation.csv");
    out_4.open("P_simulation.csv");

    for (int i=0; i<N; i++) {
        Reset(0);
        for (int j=1; j<=L; j++) {
            double S, z;
            z = rnd.Gauss(0,1);
            S = S_0 * exp((r - pow(sigma,2)/2)*T + sigma*z*sqrt(T));
            if (S > K) {
                C_blk_ave += exp(-r*T)*(S - K);
            }
            else{
                P_blk_ave += exp(-r*T)*(K-S);
            }
        }
        C_sum_prog += C_blk_ave/L;
        C_sum_prog_2 += pow(C_blk_ave/L,2);
        P_sum_prog += P_blk_ave/L;
        P_sum_prog_2 += pow(P_blk_ave/L,2);
        if(i==0){
            out_1 << C_sum_prog/(i+1) << ", " << 0 << endl;
            out_2 << P_sum_prog/(i+1) << ", " << 0 << endl;
        }
        else{
            out_1 << C_sum_prog/(i+1) << ", " << sqrt(abs(C_sum_prog_2/(i+1) - pow(C_sum_prog/(i+1),2))/i) << endl;
            out_2 << P_sum_prog/(i+1) << ", " << sqrt(abs(P_sum_prog_2/(i+1) - pow(P_sum_prog/(i+1),2))/i) << endl;
        }
    }
    Reset(1);

//  seconda parte: simulo il GBM e rifaccio
    for (int i=0; i<N; i++) {
        Reset(0);
        for (int j=1; j<=L; j++) {
            double S, z;
            z = rnd.Gauss(0,1);
            S = GBM(100, rnd, r, sigma);
            if (S > K) {
                C_blk_ave += exp(-r*T)*(S - K);
            }
            else{
                P_blk_ave += exp(-r*T)*(K-S);
            }
        }
        C_sum_prog += C_blk_ave/L;
        C_sum_prog_2 += pow(C_blk_ave/L,2);
        P_sum_prog += P_blk_ave/L;
        P_sum_prog_2 += pow(P_blk_ave/L,2);
        if(i==0){
            out_3 << C_sum_prog/(i+1) << ", " << 0 << endl;
            out_4 << P_sum_prog/(i+1) << ", " << 0 << endl;
        }
        else{
            out_3 << C_sum_prog/(i+1) << ", " << sqrt(abs(C_sum_prog_2/(i+1) - pow(C_sum_prog/(i+1),2))/i) << endl;
            out_4 << P_sum_prog/(i+1) << ", " << sqrt(abs(P_sum_prog_2/(i+1) - pow(P_sum_prog/(i+1),2))/i) << endl;
        }
    }
    out_1.close();
    out_2.close();
    out_3.close();
    out_4.close();
    return 0;
}

void Inizialization(){
    
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
}

void Reset(int a){
    if (a==0) {
        C_blk_ave=0;
        P_blk_ave=0;
    }
    else{
        C_sum_prog=0;
        C_sum_prog_2=0;
        C_blk_ave=0;
        P_sum_prog=0;
        P_sum_prog_2=0;
        P_blk_ave=0;
    }
}

double GBM(int n_steps, Random rnd, double mu, double sigma){
    double GBM = S_0;
    ofstream out;
    out.open("gbm.csv");
    out << GBM << endl;
    for (int i=0; i<n_steps; i++) {
        double z = rnd.Gauss(0,1);
        GBM = GBM * exp((mu - 0.5*pow(sigma,2))/n_steps + sigma*z/sqrt(n_steps));
        out << GBM << endl;
    }
    out.close();
    return GBM;
}
