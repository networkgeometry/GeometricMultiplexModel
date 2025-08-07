#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include "include_me/S1_realization.h"
#include <sstream>
#include <random>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>

// Custom library for specific Gaussian hypergeometric functions.
#include "include_else/hyp2f1.hpp"

using namespace std;

#include <cmath>

struct Params {
    int N = -1;
    std::string parameter_file;
    std::string correlation_file;
};

bool parse_args(int argc, char* argv[], Params &params) {
    // Flag-based parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-N" && i + 1 < argc) {
            params.N = std::stoi(argv[++i]);
        } else if (arg == "-pf" && i + 1 < argc) {
            params.parameter_file = argv[++i];
        } else if (arg == "-cf" && i + 1 < argc) {
            params.correlation_file = argv[++i];
        }
    }

    // Positional fallback if needed
    if ((params.N == -1 || params.parameter_file.empty() || params.correlation_file.empty()) and argc == 4) {
        params.N = std::stoi(argv[1]);
        params.parameter_file = argv[2];
        params.correlation_file = argv[3];
    }
    
    // Final validation
    if (params.N == -1 || params.parameter_file.empty() || params.correlation_file.empty()) {
        return false;
    }

    return true;
}

void load_parameters(string filename, vector<double>& gammas, vector<double>& betas, vector<double>& kavgs, vector<int>& seeds)
{
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(infile, line)) {
        stringstream ss(line);
        double gamma, beta, kavg;
        int seed;
        ss >> gamma >> beta >> kavg >> seed;
        gammas.push_back(gamma);
        betas.push_back(beta);
        kavgs.push_back(kavg);
        seeds.push_back(seed);
    }
    infile.close();
}

void load_correlations(string filename, vector<double>& nus, vector<double>& gs)
{
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(infile, line)) {
        stringstream ss(line);
        double nu, g;
        ss >> nu >> g;
        nus.push_back(nu);
        gs.push_back(g);
    }
    infile.close();
}

int main(int argc, char* argv[])
{   
    // Parse command line arguments
    Params params;
    if (!parse_args(argc, argv, params)) {
        std::cerr << "Usage:\n"
                  << "  ./GENNET -N 1000 -pf param.txt -cf corr.txt\n"
                  << "  or\n"
                  << "  ./GENNET 1000 param.txt corr.txt\n";
        return 1;
    }


    // Load the parameters and correlations
    std::vector<double> gammas;
    std::vector<double> betas;
    std::vector<double> kavgs;
    std::vector<double> nus;
    std::vector<double> gs;
    std::vector<int> seeds;
    load_parameters(params.parameter_file, gammas, betas, kavgs, seeds);
    load_correlations(params.correlation_file, nus, gs);

    // Initialize multiplex
    std::vector<S1_realization> multiplex;
    // Generate the first layer
    multiplex.push_back(S1_realization(gammas.at(0), betas.at(0), kavgs.at(0), params.N, seeds.at(0)));
    // Generate the subsequent layers
    for (int i = 1; i < gammas.size(); ++i) {
        multiplex.push_back(S1_realization(multiplex.at(i-1), gammas.at(i), betas.at(i), kavgs.at(i), nus.at(i-1), gs.at(i-1), seeds.at(i)));
    }

    // Output the edge files and coordinates for each layer
    std::stringstream ss;
    std::ofstream coordfile;
    for (int i = 0; i < multiplex.size(); ++i) {
        ss << "layer" << i << ".coords";
        coordfile.open(ss.str());
        ss.str(""); 
        for (int j = 0; j < params.N; ++j) {
            coordfile << multiplex.at(i).thetas.at(j) << " " << multiplex.at(i).kappas.at(j) << endl;
        }
        coordfile.close(); 
        multiplex.at(i).output_edgefile("layer" + std::to_string(i) + ".edge");
    }

    


    

  

    

    
    


}


