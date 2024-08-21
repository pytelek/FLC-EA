using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <random>
#include <chrono>

/* Ponizsze parametry zmieniac zaleznie od potrzeb */
//liczebnosc populacji
#define POPSIZE 25
//maksymalna liczba pokolen
#define MAXGENS 250
//liczba zmiennych w zadaniu
#define NVARS 50
//prawdopodobienstwo krzyzowania
#define PXOVER 0.8
//prawdopodobienstwo mutacji
#define PMUTATION 0.1
//rzeczywista liczba genow
#define NBIN 48
//delta rozkladu normalnego
#define delta 10
// liczba PI
#define M_PI 3.14159
//dlugosc historii
#define HISTORIA 5
//liczba regul
#define MAX_REGUL 8

//wartosc maksymalna (stala C)
#define MAX 40000000000

#define STDDEVMIN 0.15
#define STDDEVMAX 0.25

//mean
auto mean=0;
//stddev
auto stddev = 0.2;

int liczba_mutacji_pozytywnych[HISTORIA];
int liczba_mut_poz;
double historia_fitness[HISTORIA];
double avg_fitness;
double suma_fitness;
double suma_fit;
double srednia_fit;
double suma_mut;
// double srednia_mut;

#define TRUE 1
#define FALSE 0

//genotyp, czlonek populacji
 struct genotype
 {
  double gene[NVARS];			    //lancuch zmiennych
  double fitness;					//dopasowanie genotypu
  double upper[NVARS]; 		        //gorne ograniczenie na zmienne genotypu
  double lower[NVARS];		        //dolne ograniczenie na zmienne genotypu
  double rfitness;				    //dopasowanie wzgledne
  double cfitness;				    //dopasowanie laczne
  double jakosc;                    //jakosc w stosunku do populacji
  double wspolczynnik;              //wynik wnioskowania rozmytego
  double oblicz_dopasowanie();      //metoda obliczajaca dopasowanie
 };

//numer biezacego pokolenia
int generation;
//najlepszy osobnik
int cur_best;
//plik wyjsciowy
FILE *galog;

//populacja
genotype population[POPSIZE+1];
//nowa populacja zamieniajaca stare pokolenie
genotype newpopulation[POPSIZE+1];
//srednia funkcja przystosowania poprzedniego pokolenia
double	avg_pop;
//licznik zegara
double zegar;
//ograniczenia dolne i gorne zakresu zmiennych
double lbound, ubound;

/**************************************************************************/
// zmienne dla fuzzy logic
/**************************************************************************/
double lmp5_procent;
double avg_fitness_ratio;
double dystans;

double mf_liczba_mutacji[2];
double mf_wsp_przystosowania[2];
double mf_dystans[2];
double mf_jakosc[2];
double mf_wynik[2];

double f_mut_m, f_mut_d;
double f_przyst_h_m, f_przyst_h_d;
double f_dyst_m, f_dyst_d;
double f_jakosc[POPSIZE][2];

double r[MAX_REGUL];

double wynik_ostry;

/**************************************************************************/
// funkcja optymalizowana
/**************************************************************************/

double genotype::oblicz_dopasowanie()
{
    int k;
    double sum;
	double dopasowanie = 0;

    sum=0;

    for (k=0; k<NVARS; k++)
    {
        sum = sum + pow(1000000, ((k)/(NVARS-1))) * gene[k] * gene[k];
    }

    dopasowanie = sum;

return (MAX - dopasowanie);
};

/**************************************************************************/
//Wymaga C++11
//Settings/Compiler/Have g++ follow C++11
/**************************************************************************/
double losuj(double m, double s)
{
    double wynik;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator(seed);
    normal_distribution<double> distribution(m, s);

    wynik = distribution(generator);

    return wynik;
}

/*************************************************************************/
/*Procedura inicjujaca, nadaje genom wartosci wewnatrz ograniczen.       */
/*Ustawia na zero wartosci dopasowan dla wszystkich czlonkow popu-      */
/*lacji. Wczytuje ograniczenia dolne i gorne zmiennych z pliku           */
//wejsciowego "gadata.txt", nastepnie generuje losowo wartosci           */
//wewnatrz tych ograniczen dla kazdego genu w kazdym genotypie w         */
//populacji. Format pliku wejsciowego "gadata.txt" jak ponizej:          */
/*ograniczenie_dolne_zmiennej1  ograniczenie_gorne_zmiennej1             */
/*ograniczenie_dolne_zmiennej2  ograniczenie_gorne_zmiennej2 itd...      */
/*************************************************************************/

void initialize(void)
{
 FILE* infile;
 FILE* popfile;
 int i, j;
 double tmp;

 if ((infile=fopen("gadata.txt","r"))==NULL)
    {
     fprintf(galog,"\nNie moge otworzyc pliku wejsciowego!/n");
     exit(1);
    }
 if ((popfile=fopen("popdata.txt","r"))==NULL)
    {
     fprintf(galog,"\nNie moge otworzyc pliku popdata!/n");
     exit(2);
    }

 for (i=0; i<NVARS; i++)
    {
     fscanf(infile, "%lf", &lbound);
     fscanf(infile, "%lf", &ubound);

     for (j=0; j<POPSIZE; j++)
        {
         population[j].fitness=0;
         population[j].rfitness=0;
         population[j].cfitness=0;
         population[j].lower[i]=lbound;
         population[j].upper[i]=ubound;
        }
    }


for (i=0; i<POPSIZE; i++)
{
  for (j=0; j<NVARS; j++)
  {
  fscanf(popfile, "%lf", &tmp);
  population[i].gene[j]=tmp;
  }
}
 fclose(infile);
 fclose(popfile);


for (j=0; j<POPSIZE; j++)
  {
  population[j].wspolczynnik = 1;
  }

}

/*************************************************************************/
/*Generator liczb losowych generujacy wartosci wewnatrz                  */
/*ograniczen                                                             */
/*************************************************************************/

 double randval(double low, double high)
 {
  double val;
  val=((double)(rand()%1000/1000.0)*(high - low)+low);
  return (val);
 }

/*************************************************************************/
/*Funkcja oceny ustalana przez uzytkownika.                              */
/*Po kazdej zmianie funkcji nalezy ponownie skompilowac program.         */
/*Obecnie jest to funkcja: x[1]^2-x[1]*x[2]+x[3]                         */
/*************************************************************************/

void evaluate (void)
{
 int mem;

 for (mem=0; mem<POPSIZE; mem++)
    {
     population[mem].fitness = population[mem].oblicz_dopasowanie();

     if (population[mem].fitness < 0) population[mem].fitness = 0;
    }
}
/*************************************************************************/
/*Procedura keep_the_best. Zapamietuje ona najlepszego dotychcza-        */
/*wego czlonka populacji. Uwaga: w ostatnim elemencie macierzy           */
/*population znajduje sie kopia najlepszego osobnika.                    */
/*************************************************************************/

void keep_the_best (void)
{
 int mem;
 int i;
 cur_best =0;	//zapamietanie wskaznika najlepszego osobnika

 for (mem=0; mem<POPSIZE; mem++)
    {
     if (population[mem].fitness>population[POPSIZE].fitness)
       {
        cur_best=mem;
        population[POPSIZE].fitness=population[mem].fitness;
       }
    }

 for( i=0; i<NVARS; i++)
    {
     population[POPSIZE].gene[i]=population[cur_best].gene[i];
    }
}

/*************************************************************************/
/*Procedura elitist. Najlepszy osobnik z poprzedniego pokolenia          */
/*jest zapamietywany jako ostatni w macierzy. Jezeli najlepszy           */
/*osobnik z biezacego pokolenia jest gorszy niz najlepszy osobnik        */
/*z poprzednich pokolen, to ten ostatni zastepuje najgorszego            */
/*osobnika biezacej populacji                                            */
/*************************************************************************/

void elitist (void)
{
 int i;
 double best, worst;						//najlepsza i najgorsza wartosc dopasowania
 int best_mem=0, worst_mem=0;		    	//wskazniki do najlepszego i najgorszego osobnika

best = population[0].fitness;
 worst = population[0].fitness;
  for (i=1; i<POPSIZE; i++)
    {
    if (population[i].fitness >= best)
       {
           best = population[i].fitness;
           best_mem = i;
       }
    if (population[i].fitness <= worst)
       {
           worst=population[i].fitness;
           worst_mem=i;
       }
    }

/*************************************************************************/
/*Jezeli najlepszy osobnik z nowej populacji jest lepszy niz             */
/*najlepszy osobnik z poprzednich populacji, to skopiuj                  */
/*najlepszego z nowej populacji, jezeli nie, to zastap                   */
/*najgorszego osobnika z biezacej populacji przez najlepszego            */
/*z poprzednich pokolen                                                  */
/*************************************************************************/

     if (best>=population[POPSIZE].fitness)
       {
        for (i=0; i<NVARS; i++)
           {
            population[POPSIZE].gene[i]=population[best_mem].gene[i];
           }
        population[POPSIZE].fitness=population[best_mem].fitness;
       }
     else
       {
        for (i=0; i<NVARS; i++)
           {
            population[worst_mem].gene[i]=population[POPSIZE].gene[i];
           }
        population[worst_mem].fitness=population[POPSIZE].fitness;
       }
}

/*************************************************************************/
/*Procedura wyboru. Selekcja standardowa proporcjonalna do zadania       */
/*maksymalizacji realizujaca model elitarny - zapewnia, ze zawsze        */
/*przezywa najlepszy osobnik populacji.                                  */
/*************************************************************************/

void select (void)
{
 int mem, i, j;
 double p;
 double suma_fitness_k;
 double avg_fitness_k;

    suma_fitness=population[0].fitness;
    for (i=1; i<POPSIZE; i++)
    {
        suma_fitness += population[i].fitness;
    }
    avg_fitness = suma_fitness/POPSIZE;


    for (i=HISTORIA-2; i>=0; i--)
    {
        historia_fitness[i+1] = historia_fitness[i];
    }
    historia_fitness[0] = avg_fitness;

    suma_fitness_k=population[0].fitness * population[0].wspolczynnik;

    for (mem=1; mem<POPSIZE; mem++)
    {
        suma_fitness_k += population[mem].fitness * population[mem].wspolczynnik;
    }
    avg_fitness_k = suma_fitness_k/POPSIZE;

for (mem=0; mem<POPSIZE; mem++)
    {
     population[mem].rfitness=population[mem].fitness* population[0].wspolczynnik/suma_fitness_k;
    }

population[0].cfitness=population[0].rfitness;
 for (mem=1; mem<POPSIZE; mem++)
    {
     population[mem].cfitness=population[mem-1].cfitness+
                                population[mem].rfitness;
    }

 for( i=0; i<POPSIZE; i++)
    {
     p=rand()%1000/1000.0;
     if (p<population[0].cfitness)
       newpopulation[i]=population[0];
     else
       {
        for (j=0; j<POPSIZE; j++)
           {
            if (p>=population[j].cfitness&&p<population[j+1].cfitness)
              newpopulation[i]=population[j+1];
           }
       }
    }

 for (i=0; i<POPSIZE; i++)
    {
     population[i]=newpopulation[i];
    }
}

/*************************************************************************/
/*Funkcja pomocnicza wymieniajaca wzajem dwie zmienne.                   */
/*************************************************************************/

void swap (double *x, double *y)
{
 double temp;
 temp=*x;
 *x=*y;
 *y=temp;
}

/*************************************************************************/
/*Krzyzowanie. Wykonanie krzyzowania dwojga wybranych rodzicow.          */
/*************************************************************************/

void Xover (int one, int two)
{
 int i;
 int point;	//punkt krzyzowania

 if (NVARS>1)
   {
    if (NVARS == 2)
      {
       point = 1;
      }
    else
      {
       point = (rand() % (NVARS - 1) +1);
      }

 for (i=0; i<point; i++)
    {
     swap(&population[one].gene[i], &population[two].gene[i]);
    }
  }
}


/*************************************************************************/
/*Wybor do krzyzowania. Wybor dwojga rodzicow, ktorzy wezma udzial       */
/*w krzyzowaniu. Funkcja realizuje krzyzowanie jednopunktowe.            */
/*************************************************************************/

void crossover (void)
{
 int mem, one;
 int first=0;	//obliczenie liczby wybranych osobnikow.
 double x;
 for (mem=0; mem<POPSIZE; ++mem)
    {
     x=rand()%1000/1000.0;
     if (x<PXOVER)
       {
        ++first;
        if (first % 2 == 0)
          {
           Xover (one, mem);
          }
        else
          {
           one=mem;
          }
       }
    }
 }

/*************************************************************************/
/*Mutacja. Losowa mutacja jednorodna. Zmienna wybrana do mutacji jest    */
/*zamieniana na wartosc losowa zawarta wewnatrz ograniczenia dolnego i   */
/*gornego tej zmiennej.                                                  */
/*************************************************************************/

void mutate (void)
{

 int i, j;
 double x;

liczba_mut_poz=0;

 for (i=0; i<POPSIZE; i++)
    {
     for (j=0; j<NVARS; j++)
        {
         x=rand() % 1000/1000.0;
         if (x<PMUTATION)
           {
           population[i].gene[j]=population[i].gene[j] + losuj(mean, stddev);

           if (population[i].gene[j]<population[i].lower[j])
              {
              population[i].gene[j]=population[i].lower[j];
              }
           if (population[i].gene[j]>population[i].upper[j])
              {
              population[i].gene[j]=population[i].upper[j];
              }

           }
        }
        if (population[i].oblicz_dopasowanie() < population[i].fitness)
            {
                liczba_mut_poz++;
            }
    }
    for (i=HISTORIA-2; i>=0; i--)
    {
        liczba_mutacji_pozytywnych[i+1] = liczba_mutacji_pozytywnych[i];
    }
    liczba_mutacji_pozytywnych[0] = liczba_mut_poz;
}

/**************************************************************************/
void mf()
{
    mf_liczba_mutacji[0] = 0.25;
    mf_liczba_mutacji[1] = 0.375;

    mf_wsp_przystosowania[0] = 1.0;
    mf_wsp_przystosowania[1] = 1.1;

    mf_dystans[0] = 0.85;
    mf_dystans[1] = 0.95;

    mf_jakosc[0] = 0.9;
    mf_jakosc[1] = 1.1;

    mf_wynik[0] = 0.8;
    mf_wynik[1] = 1.2;


    //cout << "mf wczytane" << endl;
}

/**************************************************************************/
void rozmyj(double mut, double przyst_h, double dyst)
{
    int mem;

    if (mut<mf_liczba_mutacji[0])
    {
        f_mut_m=1;
        f_mut_d=0;
    }
    else
    {
        if (mut<mf_liczba_mutacji[1])
        {
            f_mut_d=(mut-mf_liczba_mutacji[0])/(mf_liczba_mutacji[1]-mf_liczba_mutacji[0]);
            f_mut_m=1-f_mut_d;
        }
        else
        {
            f_mut_m=0;
            f_mut_d=1;
        }
    }

    if (przyst_h<mf_wsp_przystosowania[0])
    {
        f_przyst_h_m=1;
        f_przyst_h_d=0;
    }
    else
    {
        if (przyst_h<mf_wsp_przystosowania[1])
        {
            f_przyst_h_d=(przyst_h-mf_wsp_przystosowania[0])/(mf_wsp_przystosowania[1]-mf_wsp_przystosowania[0]);
            f_przyst_h_m=1-f_przyst_h_d;
        }
        else
        {
            f_przyst_h_m=0;
            f_przyst_h_d=1;
        }
    }

     if (dyst<mf_dystans[0])
    {
        f_dyst_m=1;
        f_dyst_d=0;
    }
    else
    {
        if (dyst<mf_dystans[1])
        {
            f_dyst_d=(dyst-mf_dystans[0])/(mf_dystans[1]-mf_dystans[0]);
            f_dyst_m=1-f_dyst_d;
        }
        else
        {
            f_dyst_m=0;
            f_dyst_d=1;
        }
    }

 for (mem=0; mem<POPSIZE;mem++)
{
    if (population[mem].jakosc<mf_jakosc[0])
    {
        f_jakosc[mem][0]=1;
        f_jakosc[mem][1]=0;
    }
    else
    {
        if (population[mem].jakosc<mf_jakosc[1])
        {
            f_jakosc[mem][1]=(population[mem].jakosc-mf_jakosc[0])/(mf_jakosc[1]-mf_jakosc[0]);
            f_jakosc[mem][0]=1-f_jakosc[mem][1];
        }
        else
        {
            f_jakosc[mem][0]=0;
            f_jakosc[mem][1]=1;
        }
    }
}
}

/**************************************************************************/
void wyostrz()
{
    int mem;
    double m, d, wynik;

for (mem=0; mem<POPSIZE; mem++)
{
    m = f_mut_m + f_przyst_h_m + f_dyst_m + f_jakosc[mem][0];
    d = f_mut_d + f_przyst_h_d + f_dyst_d + f_jakosc[mem][1];
    wynik = d/(m+d);

    wynik_ostry = mf_wynik[0] + wynik * (mf_wynik[1]-mf_wynik[0]);

    population[mem].wspolczynnik = wynik_ostry;
}
}

/*************************************************************************/
/*Program glowny. W kazdym pokoleniu nastepuje wybor najlepszych         */
/*osobnikow, krzyzowanie i mutacja, a nastepnie ocena kolejnej populacji */
/*az do spelnienia warunku koncowego.                                    */
/*************************************************************************/

int main(int argc, char *argv[])
{
 int i;

  srand(time(0));

 if ((galog=fopen("galog.txt","a"))==NULL)
   {
    exit(3);
   }
 mf();

 generation = 0;

 initialize();

 evaluate();

 keep_the_best();

 while (population[POPSIZE].fitness<= MAX - 1)
      {
       generation++;
/**************************************************************************/
//Analiza historii fitness
/**************************************************************************/
if (generation>=HISTORIA)
{
    suma_mut=0;
    suma_fit=0;

    for (i=0; i<HISTORIA; i++)
    {
        suma_fit += historia_fitness[i];
    }
    srednia_fit = suma_fit / HISTORIA;
    avg_fitness_ratio = avg_fitness / srednia_fit;

    for (i=0; i<HISTORIA; i++)
    {
        suma_mut += liczba_mutacji_pozytywnych[i];
    }
    lmp5_procent = suma_mut / HISTORIA / POPSIZE;

    dystans = avg_fitness / population[POPSIZE].fitness;

    rozmyj(lmp5_procent, avg_fitness_ratio, dystans);
    wyostrz();
}

       select();
       crossover();
       mutate();
       evaluate();
       elitist();
    }

cout << "koniec = " << generation << endl;

fclose (galog);

 for (i=0; i<NVARS; i++)
    {
     cout << population[POPSIZE].gene[i] << " ";
    }
 cout << endl;
 cout << population[POPSIZE].fitness << endl;
 zegar=clock()/1000000.0;
 cout << "Czas wykonania [s] " << zegar << endl;


 printf ("Sukces\n");

return EXIT_SUCCESS;
}




