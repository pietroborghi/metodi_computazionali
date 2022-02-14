#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

struct cartesio {
  double x;
  double y;
  double z;
};

double dist (cartesio r1, cartesio r2, double l);

using namespace std;

int main () {

  // PARAMETRI DEL PROBLEMA
  int N;    // numero di atomi
  int N_samples, N_samples_max;  // numero di sample da analizzare
  double l;   // lunghezza box
  cartesio *atoms;  // array posizioni atomi

  int M;  // numero di punti su cui valutare g(r)
  double r_min, r_max;  // r_min e r_max su cui valutare g(r)
  double sigma;  // parametro broadening per delta di Dirac
  std::string filename_last, filename_samples, filename_out;
  std::ifstream fin_last, fin_xyz;
  char c_multisample;
  bool multisample = false;

  std::cout << "Nome file sample finale: \n";
  std::cin >> filename_last;

  std::cout << "Vuoi analizzare piu' sample? (t/f) \n";
  std::cin >> c_multisample;
  if (c_multisample == 't') {
    multisample = true;
    std::cout << "Nome file xyz samples: \n";
    std::cin >> filename_samples;
    fin_xyz.open ("./data/" + filename_samples, std::ios::in);
    fin_xyz >> N_samples_max;
    std::cout << "Numero di sample presenti: \t" << N_samples_max << '\n'
              << "Quanti sample vuoi analizzare? (oltre a LAST) \n";
    std::cin >> N_samples;
    if (N_samples > N_samples_max)
      N_samples = N_samples_max;
  }

  std::cout << "Inserire numero di punti su cui valutare pair correlation: \n";
  std::cin >> M;
  std::cout << "Inserire r_min e r_max da valutare: \n";
  std::cin >> r_min >> r_max;
  std::cout << "Inserire parametro sigma (larghezza delta Dirac): \n";
  std::cin >> sigma;
  std::cout << "Nome file output: \n";
  std::cin >> filename_out;

  // LETTURA LAST SAMPLE DA FILE
  fin_last.open ("./data/" + filename_last, std::ios::in);
  fin_last >> N >> l;
  cout << "Numero atomi: \t" << N << '\n'
       << "Lato cella simulazione: \t" << l << '\n';
  // allocazione e riempimento array dati
  atoms = new cartesio [N];
  for (int i=0; i<N; ++i)
    fin_last >> atoms[i].x >> atoms[i].y >> atoms[i].z;

  // PREPARAZIONE FILE OUTPUT
  std::ofstream fout;
  fout.open("./data/" + filename_out, std::ios::out);

  // CALCOLO PAIR CORRELATION FUNCTION
  double  r_ij;
  double *r, *g_r;        // array con r e g_r per calcolare medie
  // allocazione memoria
  r = new double [M];
  g_r = new double [M];

  for (int k=0; k<M; ++k) {

    // assegnazione valore di r
    r[k] = r_min + (r_max - r_min)*k/(M-1);
    // calcolo g_r vero e proprio
    g_r[k] = 0.0;

    for (int i=0; i<N; ++i)
      for (int j=0; j<N; ++j)
        if ( i!=j ) {
          // calcolo distanza r_ij
          r_ij = dist (atoms[i], atoms[j], l);
          // calcolo g_r (approssimando delta di Dirac a gaussiana)
          g_r[k] += 1/(sigma*sqrt(2*M_PI)) * exp ( -pow((r[k] - r_ij), 2)/(2*pow(sigma, 2)) );
        }
    g_r[k] *= pow(l, 3) / (4*M_PI*pow(r[k], 2)*pow(N, 2));

  }

  // CICLO SUI SAMPLE (se multisample)
  if (multisample) {

    // skip dati indesiderati
    double skip;
    int N_skip = (N_samples_max - N_samples)*3*N + 1;  // +1 perchè secondo numero del file è ancora N
    for (int i=0; i<N_skip; ++i)
      fin_xyz >> skip;

    // ciclo sui sample da analizzare
    for (int i=0; i<N_samples; ++i) {
      // lettura da file
      for (int j=0; j<N; ++j)
        fin_xyz >> atoms[j].x >> atoms[j].y >> atoms[j].z;
      /* cout << atoms[0].x << '\t' << atoms[0].y << '\t' << atoms[0].z << '\n'
           << atoms[N-1].x << '\t' << atoms[N-1].y << '\t' << atoms[N-1].z << '\n'; */

      // calcolo pair correlation
      for (int k=0; k<M; ++k) {
        double g_r_calc = 0.0;

        for (int i=0; i<N; ++i)
          for (int j=0; j<N; ++j)
            if ( i!=j ) {
              r_ij = dist (atoms[i], atoms[j], l);
              g_r_calc += 1/(sigma*sqrt(2*M_PI)) * exp ( -pow((r[k] - r_ij), 2)/(2*pow(sigma, 2)) );
            }
        g_r_calc *= pow(l, 3) / (4*M_PI*pow(r[k], 2)*pow(N, 2));
        g_r[k] += g_r_calc;
      }

    }

    // calcolo media
    for (int k=0; k<M; ++k)
      g_r[k] /= (N_samples+1);    // +1 è per tenere conto anche del last sample

  }

  // SCRITTURA SU FILE
  for (int k=0; k<M; ++k)
    fout << r[k] << '\t' << g_r[k] << '\n';

  // DEALLOCAZIONE MEMORIA
  delete [] atoms;
  fout.close();

  return 0;
}

// FUNZIONE CHE CALCOLA DISTANZA TRA DUE PUNTI (TENENDO CONTO DELLE REPLICHE PERIODICHE)
double dist (cartesio r1, cartesio r2, double l) {
  double dist, dist_min = 2*l;  // valore a caso, sicuramente > della vera distanza

  // ciclo su repliche periodiche
  for (int kx=-1; kx<=1; ++kx)
    for (int ky=-1; ky<=1; ++ky)
      for (int kz=-1; kz<=1; ++kz) {
        dist = sqrt ( pow(r1.x - r2.x + l*kx, 2) + pow(r1.y - r2.y + l*ky, 2) +
          pow(r1.z - r2.z + l*kz,2) );
        if (dist < dist_min)
          dist_min = dist;
      }

  return dist_min;
}
