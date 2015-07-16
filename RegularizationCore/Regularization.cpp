//
// Created by yutong pang on 7/15/15.
//

#include "Regularization.h"
#include <iostream>					// for cout etc.
#include <algorithm>				// for sort algorithm
#include <time.h>					// for random seed
#include <math.h>					// for abs()
using namespace std;

bool fitness_sort(ga_struct x, ga_struct y)
{ return (x.fitness < y.fitness); }

inline void sort_by_fitness(ga_vector &population)
{ sort(population.begin(), population.end(), fitness_sort); }

inline void swap(ga_vector *&population,
                 ga_vector *&buffer)
{
    ga_vector *temp = population;
    population = buffer;
    buffer = temp;
}



GeneticAlgorithm::GeneticAlgorithm() {
    fexact.assign(FSIZE, 0.0);
    for (int i =0; i< FSIZE; i++){
        if (i<8){
            fexact[i] = 2;
        }
        else{
            fexact[i] = 1;
        }
    }

    A.assign(FSIZE, vector<double>(FSIZE, 0.0));
    for (int i = 0; i < FSIZE; i++){
        for (int j = 0; j< FSIZE; j++){
            double s = 1.0/FSIZE * (double(j) + 0.5);
            double t = 1.0/FSIZE * (double(i) + 0.5);
            A[i][j] = 1.0/FSIZE * d/pow(pow(d, 2.0) + pow(s-t, 2.0), 1.5);
        }
    }
    b.assign(FSIZE,0.0);
    for (int i=0; i < FSIZE; i++){
        double result = 0.0;
        for (int k=0; k < FSIZE; k++){
            result += A[k][i] * fexact[k];
        }
        b[i] = result;
    }
}

void GeneticAlgorithm::runenvolotion() {
    srand(unsigned(time(NULL)));

    ga_vector pop_alpha, pop_beta;
    ga_vector *population, *buffer;
    init_population(pop_alpha, pop_beta);
    population = &pop_alpha;
    buffer = &pop_beta;

    for (int i=0; i<GA_MAXITER; i++){
        calc_fitness(*population);
        sort_by_fitness(*population);
        print_best(*population);
        mate(*population, *buffer);
        swap(population, buffer);
    }
}

void GeneticAlgorithm::init_population(ga_vector &population, ga_vector &buffer) {
    for (int i=0;  i < GA_POPSIZE; i++) {
        ga_struct structure;
        structure.fitness = 0;
        for (int i =0; i < FSIZE; i++){
            structure.freconstruct[i] = fRand(-1.0, 10.0);
        }
        population.push_back(structure);
    }
    buffer.resize(GA_POPSIZE);
}

void GeneticAlgorithm::calc_fitness(ga_vector &population) {
    for (int i = 0; i < GA_POPSIZE; i++){
        double term1 {0.0};
        double term2{0.0};
        for (int k = 0; k < FSIZE; k++){
            term2 += pow(population[i].freconstruct[k], 2.0);
        }
        for (int s =0; s < FSIZE; s++){
            double termtemp {0.0};
            for( int k=0; k < FSIZE; k++){
                termtemp += A[k][s] * population[i].freconstruct[k];
            }
            term1 += pow(termtemp - b[s], 2.0);
        }
        population[i].fitness = term1 + pow(0.1, 2.0) *term2;
    }
}



void GeneticAlgorithm::print_best(ga_vector &gav) {
    { cout << "Best: " << " (" << gav[0].fitness << ")";
        for (int i =0; i< FSIZE; i++){
            cout << " " << gav[0].freconstruct[i];
        }
        cout << endl;
    }
}

void GeneticAlgorithm::mate(ga_vector &population, ga_vector &buffer) {
    int esize = GA_POPSIZE * GA_ELITRATE, i1, i2;
    elitism(population, buffer, esize);
    // Mate the rest
    for (int i=esize; i< GA_POPSIZE; i++){
        i1 = rand() % (GA_POPSIZE / 2);
        i2 = rand() % (GA_POPSIZE / 2);
        int spos = rand() % FSIZE;
        for (int k = 0; k < FSIZE; k++){
            if (k < spos){
                buffer[i].freconstruct[k] = population[i1].freconstruct[k];
            }
            else{
                buffer[i].freconstruct[k] = population[i2].freconstruct[k];
            }
        }
        if (rand() < GA_MUTATION) mutate(buffer[i]);
    }

}

void GeneticAlgorithm::elitism(ga_vector &population, ga_vector &buffer, int esize) {
    for (int i=0; i<esize; i++) {
        buffer[i].freconstruct = population[i].freconstruct;
        buffer[i].fitness = population[i].fitness;
    }
}


void GeneticAlgorithm::mutate(ga_struct &member) {
    for (int i = 0; i< 5; i++){
        int mp = rand() % FSIZE;
        member.freconstruct[mp] = fRand(-1, 10.0);
    }
}

double GeneticAlgorithm::fRand(double fMin, double fMax) {
    double f = (rand() % 10000)/10000.0;
    return fMin + f * (fMax - fMin);
}