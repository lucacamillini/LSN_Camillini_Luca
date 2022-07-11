#include "esercizio1.1.h"

using namespace std;

int main(){
	
	//parametri	
	int M = 1000000;
	int N = 150;
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

	//definisco i vettori che contengono i dati
    double sum_prog=0, sum_prog_2=0, sum_prog_err=0, sum_prog_err_2=0;
    
//  apertura dei file
    ofstream out_ave, out_sigma, out_chi;
    out_ave.open("avereges.csv");
    out_sigma.open("sigma.csv");
    out_chi.open("chi.csv");

//	per il calcolo della media, invece che passare da un vettore di appoggio in cui salvare 1000000 numeri casuali, ne genero uno ogni volta, faccio le somme progressive
//  e scrivo su file il risultato con la relativa varianza. Così non devo occupare ram con un vector che poi butto.
	
    for (int i=1; i<=N; i++) {
        double block_mean = 0, block_test_error=0;
        for (int j=0; j<L; j++) {
            double x = rnd.Rannyu();
            block_mean += x;
            block_test_error += pow(x - 0.5,2);
        }
        sum_prog += block_mean/L;
        sum_prog_2 += pow(block_mean/L,2);
        sum_prog_err += block_test_error/L;
        sum_prog_err_2 += pow(block_test_error/L,2);

        if(i==1){
            out_ave << sum_prog/i << ", " << 0. << endl;
            out_sigma << sum_prog_err/i << ", " << 0. << endl;
        }
        else{
            out_ave << sum_prog/i << ", " << sqrt((abs(sum_prog_2/i - pow(sum_prog/i,2)))/(i-1)) << endl;
            out_sigma << sum_prog_err/i << ", " << sqrt((abs(sum_prog_err_2/i - pow(sum_prog_err/i,2)))/(i-1)) << endl;
        }
    }

//  calcolo del chi quadro
    for (int k=0; k<100; k++) {
        int counts[100];
        double chi_quad=0;
        
    //  inizializzo il vettore di conteggio
        for (int i=0; i<100; i++) {
            counts[i]=0;
        }

    //  faccio n lanci e li catalogo nei bin
        for (int j=0; j<10000; j++) {
            double random_var = rnd.Rannyu();
            double min=0, max=1;
            while (random_var>min) {
                min = min + 1./100.;
            }
            counts[int((min-1./(2*100.))*100)] += 1;
        }
    //  calcolo il chi quadro
        for (int i=0; i<100; i++) {
            chi_quad += pow(double(counts[i])-100,2)/100.;
        }
        out_chi << chi_quad << endl;
    }
    
//  scrivo il file di media a blocchi
    out_ave.close();
    out_sigma.close();
	out_chi.close();
    rnd.SaveSeed();
    
return 0;
}


