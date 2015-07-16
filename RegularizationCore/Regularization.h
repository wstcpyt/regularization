//
// Created by yutong pang on 7/15/15.
//

#ifndef REGULARIZATION_REGULARIZATION_H
#define REGULARIZATION_REGULARIZATION_H
#include <vector>
#define GA_POPSIZE 20000
#define GA_MAXITER 10000
#define GA_ELITRATE 0.1f
#define GA_MUTATION	RAND_MAX * GA_MUTATIONRATE
#define GA_MUTATIONRATE 0.7f
#define FSIZE 20
#define GA_TARGET_SIZE 2

using namespace std;
typedef vector<double > imp_vector;
typedef vector<vector<double>> operator_matrix;
typedef vector<double> measurement;
struct ga_struct
{
    imp_vector freconstruct;
    double fitness;
    ga_struct(){
        freconstruct.assign(FSIZE, 0.0);
    }
};

typedef vector<ga_struct> ga_vector;// for brevity

class GeneticAlgorithm{
public:
    GeneticAlgorithm();
    void runenvolotion();

private:
    operator_matrix A;
    double d {0.25};
    measurement b;
    vector<double> fexact;
    void init_population(ga_vector &population,
                         ga_vector &buffer );
    void calc_fitness(ga_vector &population);
    inline void print_best(ga_vector &gav);
    void mate(ga_vector &population, ga_vector &buffer);
    void elitism(ga_vector &population,
                 ga_vector &buffer, int esize );
    void mutate(ga_struct &member);
    double fRand(double fMin, double fMax);

};
#endif //REGULARIZATION_REGULARIZATION_H
