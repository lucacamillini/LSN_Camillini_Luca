//
//  main.cpp
//  
//
//  Created by luca camillini on 28/04/22.
//

#include "main.hpp"

int main(){

    out_1.open("equilibration.dat");
    out_2.open("distribution.dat");
    out_3.open("r_data_blocking.dat");
    
    Inizialization();
    
    for (int j=0; j<500; j++) {
        Evolution(1);
        out_1 << r_0 << endl;
    }
    cout << "rate accettazione equilibrazione " << (double) i/counter << endl;
    
    Reset();
    for (int j=0; j<10000; j++) {
        Evolution(1);
        out_2 << r_0 << endl;
        write_conf();
    }
    cout << "rate accettazione " << (double) i/counter << endl;
    
//  blocking
    for (int j=0; j<N; j++) {
        r_mean=0;
        Reset();
        cout << "================ blocco " << j << "================" <<endl;
        for(int h=0; h<L; h++) {
            Evolution(1);
        }
        cout << "rate di accettazione del blocco: " << (double) i/counter << endl;
        sum_prog += r_mean/L;
        sum_prog_2 += pow(r_mean/L,2);
        if(j==0){
            out_3 << j+1 << " " << sum_prog/(j+1) << " " << 0 << endl;
        }
        else{
            out_3 << j+1 << " " << sum_prog/(j+1) << " " << sqrt((abs(sum_prog_2/(j+1) - pow(sum_prog/(j+1),2)))/j) << endl;
        }
    }

    out_1.close();
    out_2.close();
    out_3.close();
    
    return 0;
}

void Evolution(unsigned int d){
    passo();
    A_0(d);
    if (rnd.Rannyu() < A) {
        for(int j=0; j<3; j++){
            r_old[j] = r_new[j];
        }
        i++;
    }
    r_0 = sqrt(r_old[0]*r_old[0] + r_old[1]*r_old[1] + r_old[2]*r_old[2])/a;
    r_mean += r_0;
    counter++;
}

void passo(){
    for (int j=0; j<3; j++) {
        r_new[j] = r_old[j] + rnd.Rannyu(-step,step);
    }
    
}

double min(double a, double b){
    double c;
    if (a<b) {
        c=a;
    }
    else{
        c=b;
    }
    return c;
}

void A_0(unsigned int d){
    double h_new = 0, h_old = 0;
    for (int j=0; j<3; j++) {
        h_old += pow(r_old[j],2);
        h_new += pow(r_new[j],2);
    }
    h_old = sqrt(h_old);
    h_new = sqrt(h_new);
    if (d==0) {
        A = min(1, p_1(h_new)/p_1(h_old));
    }
    else{
        A = min(1, p_2(h_new, r_new[2])/p_2(h_old, r_old[2]));
    }
}

double p_1(double l){
    double p = pow(a, -3./2.) / sqrt(M_PI) * exp(-l/a);
    return pow(p,2);
}

double p_2(double h, double r_z){
    double p = pow(a, -5./2.) / 8 * sqrt(2/M_PI) * h * exp(-h/(2*a)) * r_z / h;
    return pow(p,2);
}

void write_conf(){
    ofstream out;
    out.open("conf.xyz", ios::app);
    out << r_old[0]/a << " " << r_old[1]/a << " " << r_old[2]/a << endl;
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

void Reset(){
    counter = 0;
    i = 0;
}
