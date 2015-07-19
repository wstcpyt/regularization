#include <iostream>
#include "PopulationCore/Population.h"
#include <math.h>
#define GA_POPSIZE 20000
#define PRO_SIZE 20




class Fitness{
    typedef vector<vector<double>> operator_matrix;
    typedef vector<double> measurement;
    operator_matrix A;
    double d {0.25};
    vector<double> fexact;
    measurement b;
public:
    Fitness(){
        fexact.assign(PRO_SIZE, 0.0);
        for (int i =0; i< PRO_SIZE; i++){
            if (i<8){
                fexact[i] = 2;
            }
            else{
                fexact[i] = 1;
            }
        }

        A.assign(PRO_SIZE, vector<double>(PRO_SIZE, 0.0));
        for (int i = 0; i < PRO_SIZE; i++){
            for (int j = 0; j< PRO_SIZE; j++){
                double s = 1.0/PRO_SIZE * (double(j) + 0.5);
                double t = 1.0/PRO_SIZE * (double(i) + 0.5);
                A[i][j] = 1.0/PRO_SIZE * d/pow(pow(d, 2.0) + pow(s-t, 2.0), 1.5);
            }
        }
        b.assign(PRO_SIZE,0.0);
        for (int i=0; i < PRO_SIZE; i++){
            double result = 0.0;
            for (int k=0; k < PRO_SIZE; k++){
                result += A[k][i] * fexact[k];
            }
            b[i] = result;
        }
    };
    template<typename T>
    void calc_fitness(T &population){
        for (int i = 0; i < GA_POPSIZE; i++){
            double term1 {0.0};
            double term2 {0.0};
            for (int k = 0; k < PRO_SIZE; k++){
                term2 += pow(population[i].populationproperties[k], 2.0);
            }
            for (int s =0; s < PRO_SIZE; s++){
                double termtemp {0.0};
                for( int k=0; k < PRO_SIZE; k++){
                    termtemp += A[k][s] * population[i].populationproperties[k];
                }
                term1 += pow(termtemp - b[s], 2.0);
            }
            population[i].fitness = term1 + pow(0.1, 2.0) *term2;
        }
    }


};



using namespace std;int main(){

    vector<double> populationproperties(PRO_SIZE, 0.0);
    Population<vector<double>> population(populationproperties, GA_POPSIZE);
    Fitness fitness;
    population.init(-1.0, 10.0);
    for (int i=0; i < 10000; i++){
        fitness.calc_fitness((*population.populationP));
        population.sortByFitness();
        cout << "Best: " << " (" << (*population.populationP)[0].fitness << ")";
        for (int i =0; i< PRO_SIZE; i++){
            cout << " " << (*population.populationP)[0].populationproperties[i];
        }
        cout << endl;
        population.mate(0.1, 0.7);
    }


}