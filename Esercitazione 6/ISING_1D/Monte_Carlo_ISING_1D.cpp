/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
//  Simulazione a T fissata
    Input(0); //Inizialization
    ofstream eq;
    for (int i=0; i<100; i++) { //Equilibration
        eq.open("equilibration.0", ios::app);
        Move(metro);
        Measure();
        eq << temp << " " << walker[iu] << endl;
        eq.close();
    }
    
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }
        Averages(iblk,0);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    
    int counter = 0;
    for (double T=0.5; T<2.1; T+=0.1) {
        cout << "================== SIMULAZIONE A T = " << T << " ==================" << endl;
        Input(T); //Inizialization
        temp = T;
        beta = 1./temp;
        for (int i=0; i<100; i++) { //Equilibration
            Move(metro);
            Measure();
            eq.open("equilibration.0", ios::app);
            eq << temp << " " << walker[iu] << endl;
            eq.close();
        }
        counter++;
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move(metro);
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk,1);   //Print results for current block
        }
        ConfFinal(); //Write final configuration
    }

    return 0;
}


void Input(double T)
{
  ifstream ReadInput;
    if (T==0) {
        cout << "Classic 1D Ising model             " << endl;
        cout << "Monte Carlo simulation             " << endl << endl;
        cout << "Nearest neighbour interaction      " << endl << endl;
        cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
        cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    }
    
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");
    
    ReadInput >> temp;
    if (T==0) {
        beta = 1.0/temp;
        cout << "Temperature = " << temp << endl;
    }
    else{
        temp = T;
        beta = 1.0/temp;
        cout << "Temperature = " << temp << endl;
    }

  ReadInput >> nspin;
    if (T==0) {
        cout << "Number of spins = " << nspin << endl;
    }
  ReadInput >> J;
    if (T==0) {
        cout << "Exchange interaction = " << J << endl;
    }
    ReadInput >> h;
    if (T==0) {
        cout << "External field = " << h << endl << endl;
    }
    ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
    if(T==0){
      if(metro==1) cout << "The program perform Metropolis moves" << endl;
      else cout << "The program perform Gibbs moves" << endl;
      cout << "Number of blocks = " << nblk << endl;
      cout << "Number of steps in one block = " << nstep << endl << endl;
    }
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
    if (T!=0) {
        ifstream ReadConfig;
        ReadConfig.open("config.final");
        int i=0;
        while (!ReadConfig.eof()) {
            ReadConfig >> s[i];
            i++;
        }
    }
    else{
      for (int i=0; i<nspin; ++i)
      {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      }
    }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
    if(T==0){
        cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
    }
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        energy_old = Boltzmann(s[o],o);
        energy_new = Boltzmann(-s[o],o);
        p = exp(-beta*(energy_new-energy_old));
        if (rnd.Rannyu()<min(1,p)) {
            s[o]=-s[o];
            accepted++;
        }
        attempted++;
    }
    else //Gibbs sampling
    {
        energy_up = Boltzmann(1,o);
        energy_down = Boltzmann(-1,o);
        p = exp(-beta*(energy_down-energy_up));
        if (rnd.Rannyu()<1/(1+p)) {
            s[o]=1;
            accepted++;
        }
        else{
            s[o]=-1;
            accepted++;
        }
        attempted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
      m += s[i];
  }
    walker[iu] = u;
    walker[ic] = u*u;
    walker[im] = m;
    walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk,int h) //Print results for current block
{
    
    ofstream Ene, Heat, Mag, Chi;
    if (h==0) {
        cout << "Block number " << iblk << endl;
        cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
    
    stima_u = blk_av[iu]/blk_norm/(double) nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if(h==0){
        Ene.open("output.ene.0",ios::app);
        Ene << iblk << " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
        Ene.close();
    }
    if (iblk%nblk==0 && h==1) {
        ofstream out_1;
        out_1.open("ave_output.ene.0",ios::app);
        out_1 << " " << temp << " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
        out_1.close();
    }
    
    stima_m = abs(blk_av[im]/blk_norm/(double)nspin); //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    if(h==0){
        Mag.open("output.mag.0",ios::app);
        Mag << iblk << " "  << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
        Mag.close();
    }
    if (iblk%nblk==0 && h==1) {
        ofstream out_2;
        out_2.open("ave_output.mag.0",ios::app);
        out_2 << temp << " "  << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
        out_2.close();
    }
    
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin; //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    if(h==0){
        Heat.open("output.heat.0",ios::app);
        Heat << iblk << " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
        Heat.close();
    }
    if (iblk%nblk==0 && h==1) {
        ofstream out_3;
        out_3.open("ave_output.heat.0",ios::app);
        out_3 << temp << " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
        out_3.close();
    }
    
    stima_x = blk_av[ix]/blk_norm/(double)nspin;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    if(h==0){
        Chi.open("output.chi.0",ios::app);
        Chi << iblk << " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
        Chi.close();
    }
    if (iblk%nblk==0 && h==1) {
        ofstream out_4;
        out_4.open("ave_output.chi.0",ios::app);
        out_4 << temp << " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
        out_4.close();
    }
    if (h==0) {
        cout << "----------------------------" << endl << endl;
    }
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double min(double a,double b){
    if (a<b) {
        return a;
    }
    else{
        return b;
    }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
