//
//  main.cpp
//  
//
//  Created by luca camillini on 07/07/22.
//

#include "main.hpp"

using namespace std;

int main(int argc, char* argv[]){
    
//  inizializzo MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat0, stat1;
    
    if (rank == 0) {
        cout << "=============Salesman Problem - GENETIC ALGORITHM and PARALLEL COMPUTING=============" << endl;
    }
    
//  inizializzo il generatore di numeri casuali
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for (int i=0; i<rank+1; i++) {
            Primes >> p1 >> p2 ;
        }
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
    
//  inizializzo la classe (città nel quadrato)
    read_file("American_capitals.dat");
    int N_cities = 50;
    ga gen_al(N_cities,1,200, rnd, 0.5, 0.1);
    gen_al.Set_World("cities_coordinates.input");
    cout << "settate le coordinate" << endl;
    gen_al.Sort();
    gen_al.print("first_fitness_" + to_string(rank) + ".dat","first_population_" + to_string(rank) + ".dat");
    gen_al.print_half("first_best_half_" + to_string(rank) + ".dat");
    
    int chrom_0[N_cities], chrom_1[N_cities];
    int rank_0, rank_1;
    
    for (int i=0; i<size; i++) {
        if (rank==i) {
            gen_al.print_world("city_coordinates" + to_string(rank) + ".dat");
        }
    }
    
    for(int i=0; i<1000; i++){
        gen_al.Next_Gen();
        gen_al.Sort();
        //ogni 50 generazioni faccio un mix della popolazione
        if((i+1)%50==0){
            
            int itag0 = 1;
            if(rank == 0) cout << "Generazione " << i+1 << endl;

            //scelgo due rank in maniera casuale
            if(rank==0){
                rank_0 = int(rnd.Rannyu(0,size));
                rank_1 = int(rnd.Rannyu(0,size));
                while (rank_0==rank_1) {
                    rank_1 = int(rnd.Rannyu(0,size));
                }
            }
            //Mando i rank_0 e 1 in tutti i core
            MPI_Bcast(&rank_0, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
            MPI_Bcast(&rank_1, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
            
            //scelgo i due best delle due routine
            
            if (rank==rank_0) {
                vector<int> best_0;
                best_0 = gen_al.Get_Chrom(0);
                for (int j=0; j<N_cities; j++) chrom_0[j] = best_0[j];
            }
            if (rank==rank_1) {
                vector<int> best_1;
                best_1 = gen_al.Get_Chrom(0);
                for (int j=0; j<N_cities; j++) chrom_1[j] = best_1[j];
            }
            
            //Scambio i cromosomi tra i due rank selezionati
            if (rank==rank_1) {
                MPI_Send(&chrom_1[0], N_cities, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD);
                MPI_Recv(&chrom_0[0], N_cities, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD, &stat0);
                
                //Assegno il nuovo cromosoma nella popolazione
                vector<int> best;
                for (int j=0; j<N_cities; j++) best.push_back(chrom_0[j]);
                gen_al.Set_Chrom(best,0);
            }
            if(rank==rank_0){
                MPI_Send(&chrom_0[0], N_cities, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD);
                MPI_Recv(&chrom_1[0], N_cities, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD, &stat0);
                
                vector<int> best;
                for (int j=0; j<N_cities; j++) best.push_back(chrom_1[j]);
                gen_al.Set_Chrom(best,0);
            }
            gen_al.Sort();
        }
        gen_al.print("fitness_" + to_string(rank) + ".dat", "population_" + to_string(rank) + ".dat");
        gen_al.print_half("best_half_" + to_string(rank) + ".dat");
    }
    MPI_Finalize();
    return 0;
}

void read_file(string a){
    ifstream ReadInput;
    ReadInput.open(a);
    
    ofstream out;
    out.open("cities_coordinates.input");
    int i=0;
    while ( i<50 ) {
        i++;
        string city, state;
        double longitude, latitude;
        ReadInput >> city;
        ReadInput >> state;
        ReadInput >> longitude;
        ReadInput >> latitude;
        out << longitude << " " << latitude << endl;
    }
    ReadInput.close();
    out.close();
}
