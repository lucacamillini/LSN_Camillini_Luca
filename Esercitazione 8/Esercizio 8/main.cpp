//
//  main.cpp
//  
//
//  Created by luca camillini on 12/05/22.
//

#include "main.hpp"

int main() {
    out_1.open("equilibration.dat");
    out_2.open("campionamento.dat");
    out_3.open("data.dat");
    out_4.open("annealing.dat");
    
    Initialization();
//Equilibrazione
    for (int j=0; j<10000; j++) {
        Move();
        Energy();
    }
    for (int j=0; j<10000; j++) {
        Move();
        Energy();
        out_1 << j << " " << E << endl;
        out_2 << x_old << endl;
    }
    cout << "rate di accettazione: " << (double) i/counter << endl;
    Reset();
    for (int j=0; j<N; j++) {
        blk_ave=0;
        for (int h=0; h<L; h++) {
            Move();
            Energy();
        }
        sum_prog += blk_ave/L;
        sum_prog_2 += pow(blk_ave/L,2);
        if(j==0){
            out_3 << j+1 << " " << sum_prog/(j+1) << " " << 0. << endl;
        }
        else{
            out_3 << j+1 << " " << sum_prog/(j+1) << " " << sqrt((abs(sum_prog_2/(j+1) - pow(sum_prog/(j+1),2)))/j) << endl;
        }
    }
    
    mu_old = mu;
    sigma_old = sigma;
    energy_old = sum_prog/N;
    sigma_energy_old = sqrt((abs(sum_prog_2/N - pow(sum_prog/N,2)))/(N-1));
    
//  Annealing
//  Innanzitutto campiono a diverse temperature, sempre più basse
    /*double T = T_max;
    int r = 0;
    double step_mu = 0.01, step_sigma = 0.005;
    while(T > 0.001){
        r++;
        double beta = 1./T;
        
//      Ora vario in maniera casuale i parametri della distribuzione, mu e sigma, per vedere se trovo di meglio Ripeto questa operazione un po' di volte
        Reset();
        mu = abs(mu_old + 0.01*rnd.Rannyu(-1,1));
        sigma = abs(sigma_old + 0.005*rnd.Rannyu(-1,1));
    
        eval_energy();
        
        if (energy_new < -0.45) {
            step_mu = 0.005;
            step_sigma = 0.001;
        }
        
        A = min(1, exp(-beta*(energy_new-energy_old)));
        
        if (rnd.Rannyu() < A) {
            mu_old = mu;
            sigma_old = sigma;
            energy_old = energy_new;
            sigma_energy_old = sigma_energy_new;
            i++;
        }
        counter++;
        
//      Salvo su file i parametri e l'energia trovata per il GS
        out_4 << r << " " << T << " " << mu_old << " " << sigma_old << " " << energy_old << " " << sigma_energy_old << endl;
        if (r%5==0) {
            cout << "simulazione a T = " << T << endl;
            cout << "Energia media " << energy_old << endl;
            //cout << "rate accettazione " << (double) i/counter << endl;
            cout << "---------------------------------" << endl;
        }
        T = T*0.9995;
    }
    cout << mu_old << " " << sigma_old << endl;
    //  Produco dei dati per fare un check che il campionamento sia corretto
    */
    out_1.close();
    out_2.close();
    out_3.close();
//    out_4.close();
    
    return 0;
}

void Move(){
    x = x_old + rnd.Rannyu(-step, step);
    A_0();
    if (rnd.Rannyu() < A) {
        x_old = x;
        i++;
    }
    counter++;
}

void Energy(){
    double alpha1 = pow(x_old-mu,2)/pow(sigma,2)*0.5;
    double alpha2 = pow(x_old+mu,2)/pow(sigma,2)*0.5;
    double a_1 = pow(mu,2) - pow(sigma,2) + pow(x_old,2) - 2*mu*x_old;
    double a_2 = pow(mu,2) - pow(sigma,2) + pow(x_old,2) + 2*mu*x_old;
    double v = pow(x_old,4) - 5./2.*pow(x_old,2);
    
    E = -0.5/(psi(x_old)*pow(sigma,4))*(exp(-alpha1)*a_1 + exp(-alpha2)*a_2) + v;
    
    blk_ave += E;
}

void eval_energy(){
    Reset();
    for (int j=0; j<N; j++) {
        blk_ave = 0;
        for(int k=0; k<L; k++){
            Move();
            Energy();
        }
        sum_prog += blk_ave/L;
        sum_prog_2 += pow(blk_ave/L,2);
    }
    energy_new = sum_prog/N;
    sigma_energy_new = sqrt((abs(sum_prog_2/N - pow(sum_prog/N,2)))/(N-1));
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

void A_0(){
    double p_old = pow(psi(x_old),2);
    double p_new = pow(psi(x),2);
    A = min(1, p_new/p_old);
}

void Initialization(){
    //Inizializzo il generatore di numeri casuali
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
    i = 0;
    counter = 0;
    sum_prog = 0;
    sum_prog_2 = 0;
    blk_ave = 0;
}

double psi(double x){
    double alpha_1, alpha_2;
    alpha_1 = 0.5*pow(x - mu,2)/pow(sigma, 2);
    alpha_2 = 0.5*pow(x + mu,2)/pow(sigma,2);
    return exp(-alpha_1) + exp(-alpha_2);
}
