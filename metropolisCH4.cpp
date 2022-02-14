#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <random>

double const kB = 1.38064852e-23;     // Boltzmann constant
double const A = 1e-10;               // Angstrom

struct cartesio {
  double x;
  double y;
  double z;
};


void initial_conf (int N, double l, cartesio * r0);
double pot_energy (int N, double sigma, double eps, double r_co, double l,
                   const cartesio * r0);
double change_config (int N, double sigma, double eps, double r_co, double l,
                      const cartesio * r0, cartesio * r1);

using namespace std;

// random number generator
static default_random_engine e;
static uniform_real_distribution <double> random_number(0,1);

int main (int argc, const char * argv[]) {

  // PARAMETRI MOLECOLA METANO

  double rho, r_co, T;  // densità, raggio di cut-off (per potenziale LJ) e temperatura
  double eps = 140.42; // dispersion energy (in units of kB)
  double sigma = 4.015; // size of the particle (in Angrstrom)
  double m = 16.04*1.6747e-27; //massa molecola metano

  int N_steps; // numero passi simulazione
  int N_print_ene; // numero passi scrittura su file energie
  int N_print_xyz; // numero passi scrittura su file posizioni
  int N; // numero atomi
  cartesio *r0, *r1; // array posizioni atomi
  double U_pot; // energia potenziale

  std::cout << "Inserire numero di atomi: \n";
  std::cin >> N;
  std::cout << "Inserire densita': (kg/m^3) \n";
  std::cin >> rho;
  std::cout << "Inserire temperatura: (K) \n";
  std::cin >> T;
  std::cout << "Inserire raggio di cut off: (in unita' di sigma) \n";
  std::cin >> r_co;
  std::cout << "Inserire numero di passi simulazione: \n";
  std::cin >> N_steps;

  // preparazione scrittura su file
  char c_print_en, c_print_xyz;
  bool print_en, print_xyz;
  print_en = print_xyz = false;
  std::string filename_ene, filename_xyz;
  std::ofstream fout_ene, fout_xyz, fout_last;

  std::cout << "Scrivere energie su file? (t/f) \n";
  std::cin >> c_print_en;
  if (c_print_en == 't') {
    print_en = true;
    std::cout << "Inserire nome file energie: \n";
    std::cin >> filename_ene;
    fout_ene.open ("./data/" + filename_ene, std::ios::out);
    std::cout << "Stampa ogni : \n";
    std::cin >> N_print_ene;
  }
  std::cout << "Scrivere posizioni atomi su file? (t/f)  \n";
  std::cin >> c_print_xyz;
  if (c_print_xyz == 't') {
    print_xyz = true;
    std::cout << "Inserire nome file posizioni: \n";
    std::cin >> filename_xyz;
    fout_xyz.open ("./data/" + filename_xyz, std::ios::out);
    std::cout << "Stampa ogni : \n";
    std::cin >> N_print_xyz;
    // print numero samples da scrivere su file + quanti atomi per sample
    fout_xyz << (N_steps / N_print_xyz) << '\t' << N << '\n';
  }

  // CALCOLO PARAMETRI UTILI
  r_co *= sigma; // cut off radius
  double l = pow(m*N/rho, 1./3.) / A; // lato cella simulazione (in Angstrom)
  double beta = 1/T; // beta (a meno di un kB)

  std::cout << '\n' << "Lato cella Simulazione (A): " << l << '\n';

  // ALLOCAZIONE ARRAY
  r0 = new cartesio [N];
  r1 = new cartesio [N];

  // GENERAZIONE CONFIGURAZIONE DI PARTENZA
  initial_conf (N, l, r0);
  U_pot = pot_energy (N, sigma, eps, r_co, l, r0) * kB/N;      // diviso N è per avere energia per sito
  std::cout << "Energia potenziale iniziale:" << '\t' << U_pot << '\n';

  // CICLO SUL NUMERO DI PASSI
  int N_notify = N_steps/20;
  int status;
  double r, delta_en, p_ratio;

  for (int i=0; i<N_steps; ++i) {

    // aggiornamento user su stato di avanzamento
    if ( (i % N_notify) == 0 ) {
      status = floor ((i*100.0) / N_steps);
      cout << "Stato di avanzamento Metropolis: " << status << "%" << '\n';
    }

    // ALGORITMO DI METROPOLIS

    // generazione configurazione di prova e calcolo rapporto tra probabilità
    delta_en = change_config (N, sigma, eps, r_co, l, r0, r1);
    p_ratio = exp (-beta*delta_en);

    // scelta nuova configurazione
    if (p_ratio > 1) {
      for (int k=0; k<N; ++k)
        r0[k] = r1[k];
    }
    else {
      r = random_number(e);
      if (r < p_ratio)
        for (int k=0; k<N; ++k)
          r0[k] = r1[k];
    }

    // SCRITTURA SU FILE
    // scrittura energie
    if (print_en)
      if ( (i%N_print_ene == 0) ) {
        U_pot = pot_energy (N, sigma, eps, r_co, l, r0) * kB/N;
        fout_ene << U_pot <<'\n';
      }
    // scrittura posizioni
    if (print_xyz)
      if ( (i%N_print_xyz == 0) )
        for (int j=0; j<N; ++j)
          fout_xyz << r0[j].x << '\t' << r0[j].y << '\t' << r0[j].z << '\n';
  }

  U_pot = pot_energy (N, sigma, eps, r_co, l, r0) * kB/N;
  std::cout << "Energia potenziale finale:" << '\t' << U_pot << '\n';

  // scrittura su file ultima config
  fout_last.open("./data/last.txt");
  fout_last << N << '\t' << l << '\n';
  for (int i=0; i<N; ++i)
    fout_last << r0[i].x << '\t' << r0[i].y << '\t' << r0[i].z << '\n';
  fout_last.close();

  // DEALLOCAZIONE ARRAY
  delete [] r0;
  delete [] r1;

  return 0;
}

// FUNZIONI GLOBALI UTILIZZATE:

// FUNZIONE CHE CREA CONFIGURAZIONE INIZIALE: (random)

void initial_conf (int N, double l, cartesio *r0) {
  double r;
  for (int i=0; i<N; ++i) {
    r = random_number(e);
    r0[i].x = (r*l);
    r = random_number(e);
    r0[i].y = (r*l);
    r = random_number(e);
    r0[i].z = (r*l);
  }

  return;
}

// FUNZIONE CHE CALCOLA ENERGIA POTENZIALE DI UNA CONFIGURAZIONE

double pot_energy (int N, double sigma, double eps, double r_co, double l,
                   const cartesio * r0) {

  double en = 0.0; // energia potenziale
  double dist, dist_min; // distanza e distanza minima con repliche periodiche

  // ciclo su tutti gli atomi
  for (int i=0; i<N; ++i)
    for (int j=0; j<i; ++j) {  // i<j serve a evitare fattore 0.5

      dist_min = r_co;
      // ciclo su copie periodiche
      for (int kx=-1; kx<=1; ++kx)
        for (int ky=-1; ky<=1; ++ky)
          for (int kz=-1; kz<=1; ++kz) {
            dist = sqrt ( pow(r0[i].x - r0[j].x + l*kx, 2) + pow(r0[i].y - r0[j].y + l*ky, 2) +
              pow(r0[i].z - r0[j].z + l*kz,2) );
            if (dist < dist_min)
              dist_min = dist;
          }

      if (dist_min < r_co)
        en += 4.*eps*( pow(sigma/dist_min, 12) - pow(sigma/dist_min, 6) );
    }

  return en;
}

// FUNZIONE CHE CAMBIA CONFIGURAZIONE (e calcola differenza di energia)

double change_config (int N, double sigma, double eps, double r_co, double l,
                      const cartesio * r0, cartesio * r1) {

  // copia r0 su r1
  for (int i=0; i<N; ++i)
    r1[i] = r0[i];

  double r;
  int icount;
  double dist0, dist1, distMin0, distMin1;
  double delta_en = 0.0;

  // cambio configurazione
  r = random_number(e);
  icount = floor (r*N);    // atomo di cui cambiare posizione
  r = random_number(e);
  r1[icount].x = (r*l);
  r = random_number(e);
  r1[icount].y = (r*l);
  r = random_number(e);
  r1[icount].z = (r*l);

  // calcolo differenza di energia

  for (int j=0; j<N; ++j) {  // ciclo su tutti atomi
    distMin0 = distMin1 = r_co;
    // ciclo su copie periodiche
    for (int kx=-1; kx<=1; ++kx)
      for (int ky=-1; ky<=1; ++ky)
        for (int kz=-1; kz<=1; ++kz) {
          // calcolo distanza minima nuova configurazione
          dist1 = sqrt ( pow(r1[icount].x - r1[j].x + l*kx, 2) + pow(r1[icount].y - r1[j].y + l*ky, 2) +
            pow(r1[icount].z - r1[j].z + l*kz,2) );
          if (dist1 < distMin1)
            distMin1 = dist1;
          // calcolo distanza minima vecchia configurazione
          dist0 = sqrt ( pow(r0[icount].x - r0[j].x + l*kx, 2) + pow(r0[icount].y - r0[j].y + l*ky, 2) +
            pow(r0[icount].z - r0[j].z + l*kz,2) );
          if (dist0 < distMin0)
            distMin0 = dist0;
        }

    // contributo a variazione di energia di nuova configurazione
    if ( (0 < distMin1) && (distMin1 < r_co) )
      delta_en += 4.*eps*( pow(sigma/distMin1, 12) - pow(sigma/distMin1, 6) );

    // contributo a variazione di energia di vecchia configurazione
    if ( (0 < distMin0) && (distMin0 < r_co) )
      delta_en -= 4.*eps*( pow(sigma/distMin0, 12) - pow(sigma/distMin0, 6) );
  }

  return delta_en;
}
