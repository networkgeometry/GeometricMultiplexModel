#ifndef S1_REALIZATION_INCLUDED
#define S1_REALIZATION_INCLUDED

#include <vector>
#include <random>
#include <set>

#include <fstream>
#include <algorithm>

//#include "hyp2f1.hpp"

using namespace std;

class network{
    public:
        int N=0;
        int kavg;
        vector<vector<int>> edges;
        vector<int> degree;
    public:
        void output_edgefile(string name);      
};

class S1_realization: public network{
    public:
        double beta;
        double gamma;
        double mu;
        double RH2 = 0;
        int seed = 1534124;
        double kappa_0;

        vector<long double> thetas;
        vector<double> kappas;
        vector<double> rs;

    private:
        //for pareto distr.
        double kappa_c;
        //for connecting;
        long double pij;
        double xij;
        double random;

    public:
        S1_realization(double gamma, double beta, int kavg, int N, int seed);
        S1_realization(S1_realization layer1, double gamma2, double beta2, int kavg2, double nu, double g, int seed);
        S1_realization();
        void Random_Connect();
        void Fill_kappas();
        void Calc_rs();
        void Calc_RH2();
        void Fill_thetas();
        void Calc_mu();
        void calculate_numerical_mu();
};

#endif // S1_REALIZATION_INCLUDED
