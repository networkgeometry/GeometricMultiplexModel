#define _USE_MATH_DEFINES

#include "S1_realization.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <set>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/erf.hpp>

//#include <algorithm>
//#include <random>

// Custom library for specific Gaussian hypergeometric functions.
#include "../include_else/hyp2f1.hpp"

using namespace std;

#include <cmath>

//////////////////////////////////////////////////////////////////
//                                                              //
//                                                              //
//               Reading and writing edge files                 //
//                                                              //
//                                                              //
//////////////////////////////////////////////////////////////////




/*

Outputting an edge file

*/
void network::output_edgefile(string name){
    ofstream outfile;
    outfile.open(name);
    for (int i=0; i<edges.at(0).size(); i++){
        outfile << edges.at(0).at(i) << " " << edges.at(1).at(i) << endl;
    }
    outfile.close();
}

//////////////////////////////////////////////////////////////////
//                                                              //
//                                                              //
//                  Network Analysis Functions                  //
//                                                              //
//                                                              //
//////////////////////////////////////////////////////////////////

/*

Filling the master equation

*/
void network::fill_master(){
    if (master.size()!=0){
        return;
    }

	master.clear();
	master.resize(2*edges.at(1).size());
	pi.push_back(master.data());
	pf.push_back(master.data()-1);
	//initialize pointers
	for (int i=1; i<degree.size(); i++){
		pi.push_back(pi.back()+degree.at(i-1));
		pf.push_back(pi.at(i)-1);
	}
	for (int i=0; i<edges.at(0).size();i++){
		pf.at(edges.at(0).at(i)-1)=pf.at(edges.at(0).at(i)-1)+1;
		*pf.at(edges.at(0).at(i)-1)=edges.at(1).at(i);
		pf.at(edges.at(1).at(i)-1)=pf.at(edges.at(1).at(i)-1)+1;
		*pf.at(edges.at(1).at(i)-1)=edges.at(0).at(i);
	}
}

//////////////////////////////////////////////////////////////////
//                                                              //
//                                                              //
//                      S1 functions                            //
//                                                              //
//                                                              //
//////////////////////////////////////////////////////////////////


/*

Initializing an empty S1 realization

*/
S1_realization::S1_realization(){
}
/*

Initializing a filled S1_realization

*/
S1_realization::S1_realization(double gamma, double beta, int kavg, int N, int seed){
    this->gamma=gamma;
    this->beta=beta;
    this->kavg=kavg;
    this->N=N;
    this->seed=seed;

    //
    Fill_thetas();

    if (gamma>2){
        Fill_kappas();
    }
    else{
        for (int i=0; i<N; i++){
            kappas.push_back(kavg);
        }
    }

    //Calculate mu (here defined as the
    //Calc_mu();
    calculate_numerical_mu();

    Random_Connect();

}

/*

Initializing a subsequent layer of the S1 multiplex

*/

S1_realization::S1_realization(S1_realization layer1, double gamma2, double beta2, int kavg2, double nu, double g, int seed){
    this->gamma=gamma2;
    this->beta=beta2;
    this->kavg=kavg2;
    this->seed=seed;
    N = layer1.N;

    mt19937 gen{static_cast<uint32_t>(seed)};
    uniform_real_distribution<> uni{0,1};
    double a, fnu, phi1, phi2, kappa2;

    if (gamma2>2){
        kappa_0=(1-pow(N,-1))/(1-pow(N,(2-gamma)/(gamma-1)))*(gamma-2)/(gamma-1)*kavg;
        kappa_c=kappa_0*pow(N,1/(gamma-1));
        if (nu > 0 && nu < 1){
            for (int i=0; i<N; i++){
                phi1 =  std::log(1-std::pow(layer1.kappa_0/layer1.kappa_c,layer1.gamma-1));
                phi1 -= std::log(1-std::pow(layer1.kappa_0/layer1.kappas.at(i),layer1.gamma-1));
                a = std::pow(phi1,nu / (1 - nu));
                a *= 1-std::pow(layer1.kappa_c/layer1.kappa_0,1-layer1.gamma);
                a /= 1-std::pow(layer1.kappas.at(i)/layer1.kappa_0,1-layer1.gamma);
                fnu = nu / (1-nu) * boost::math::lambert_w0((1-nu)/nu*std::pow(uni(gen) / a,(nu-1)/nu));
                phi2 = std::pow(std::pow(fnu,1.0/(1-nu)) - std::pow(phi1,1.0/(1-nu)),1-nu);
                kappas.push_back(std::pow(1-(1-std::pow(kappa_0/kappa_c,gamma-1))*std::exp(-phi2),1.0/(1-gamma))*kappa_0);
            }
        }
        else if (nu == 0){
            Fill_kappas();
        }
        else if (nu == 1){
            for (int i=0; i<N; i++){
                kappa2 = 1-std::pow(layer1.kappas.at(i)/layer1.kappa_0,1-layer1.gamma);
                kappa2 *= 1-std::pow(kappa_c/kappa_0,1-gamma);
                kappa2 /= 1-std::pow(layer1.kappa_c/layer1.kappa_0,1-layer1.gamma);
                kappa2 *= -1;
                kappa2 += 1;
                kappa2 = std::pow(kappa2,1.0/(1-gamma));
                kappa2 *= kappa_0;
                kappas.push_back(kappa2);
            }
        }
    }
    else{
        for (int i=0; i<N; i++){
            kappas.push_back(kavg2);
        }
    }

    double Phia,Phib,sigma,l,theta_new;

    if (g > 0 && g < 1){
        sigma = std::min(100.,N/(4*M_PI))*(1/g - 1);
        Phib = boost::math::erf(N / (2 * sigma));
        Phia = boost::math::erf(-N / (2 * sigma));

        for (int i = 0; i < N; i++){
            l = sigma * boost::math::erf_inv((Phib- Phia)*uni(gen)+Phia);
            theta_new = std::fmod(layer1.thetas.at(i)+2*M_PI*l/N,2*M_PI);
            if (theta_new < 0){
                theta_new += 2*M_PI;
            }
            thetas.push_back(theta_new);
        }
    }
    else if (g == 0){
        Fill_thetas();
    }
    else if (g == 1){
        for (int i = 0; i < N; i++){
            thetas.push_back(layer1.thetas.at(i));
        }
    }

    calculate_numerical_mu();
    Random_Connect();
}

/*

Functions needed to build a realization

*/

//Calculating mu
void S1_realization::Calc_mu(){
    double avg_kappa(0);
    for (int i=0; i < N; i++){
        avg_kappa += kappas.at(i)/N;
    }

    if (beta>1){
        mu=beta*sin(M_PI/beta)/(2*M_PI*kavg);
    }
    else if (beta == 1){
        mu=1./(2*avg_kappa*log(N));
    }
    else if (beta>0){
        mu=pow(N,beta-1)*(1-beta)/(pow(2,beta)*avg_kappa);
    }
    else{
        mu= 1. / (avg_kappa * N);
    }
}

//Fill kappas
void S1_realization::Fill_kappas(){
    mt19937 gen{static_cast<uint32_t>(seed)};
    uniform_real_distribution<> kappa_distrib{0,1};


    kappa_0=(1-pow(N,-1))/(1-pow(N,(2-gamma)/(gamma-1)))*(gamma-2)/(gamma-1)*kavg;
    kappa_c=kappa_0*pow(N,1/(gamma-1));
    //kappa_0=(gamma-2)/(gamma-1)*kavg;
    for (int i=0; i<N; i++){
        //kappas.push_back(kappa_0*pow(1-kappa_distrib(gen),1/(1-gamma)));
        kappas.push_back(kappa_0*pow(1-kappa_distrib(gen)*(1-pow((kappa_c/kappa_0),1-gamma)),1/(1-gamma)));
    } 
}


void S1_realization::Calc_RH2(){
    if (beta >= 1){
        RH2 = 2.0 * std::log(N / M_PI);
    }
    else{
        RH2 = 2.0 * beta * std::log(N / M_PI);
    }
    RH2 -= 2 * std::log(mu * kappa_0 * kappa_0);
}

void S1_realization::Calc_rs(){
    if (RH2 == 0) {Calc_RH2();}
    for (auto it = kappas.begin(); it != kappas.end(); ++it){
        rs.push_back(RH2 - 2 * std::log(*it / kappa_0));
    }
}

void S1_realization::calculate_numerical_mu()
{


    //Newtons method to find the value of mu where in the model the expected degree is the average degree
    std::vector<double>::iterator it1,it2;

    bool keep_going(true);
    int cnt(0);
    int k1,k2, Nk1, Nk2;
    double int1,int2,tmp;
    double a,b,b1,e,eprime,emax,e1,e2;
    double kappa_prime, curly_c;
    double xi;

    if (gamma > 2){
        // Initialize mu by its value for N -> Infinity and calculate prefactors (note that in a the emax is missing cause we are calculating here <k> not <E>)
        if (beta>=1){
            mu = std::max(beta * std::sin(M_PI/beta) / (2*M_PI*kavg),0.01);
            b = (2-gamma);
            b1 = (2-gamma);
            a = 2*pow((gamma-1)/(2-gamma),2)*pow(kappa_0,2) * 1.0 / (1 - 1.0 / N) * 1.0 / (1 - 1.0 / N);
            emax = log(N) - log(2) - 2 * log(kappa_0);
            e1 = log(N) - log(2) - log(kappa_0) - log(kappa_c);
            e2 = log(N) - log(2) - 2 * log(kappa_c);
        }
        else if (beta > 0){
            mu = (1-beta)/(std::pow(2,beta)*kavg*std::pow(N,1-beta));
            b = beta+1-beta*gamma;
            b1 = 1 + 1/beta - gamma;
            a = 2*beta*beta*pow((gamma-1)/(beta+1-beta*gamma),2)*pow(kappa_0,2.0/beta) * 1.0 / (1 - 1.0 / N) * 1.0 / (1 - 1.0 / N);
            emax = log(N) - log(2) - 2 / beta * log(kappa_0);
            e1 = log(N) - log(2) - 1 / beta * log(kappa_0) -1 / beta * log(kappa_c);
            e2 = log(N) - log(2) - 2 / beta * log(kappa_c);
        }
        else{
            mu = 1.0 / (kavg * N);
            curly_c = (gamma-1) * std::pow(kappa_0, (gamma - 1)) * 1.0 / (1 - 1.0 / N);
        }

        //cout << a << " " << b << " " << emax << " " << kappa_0 << endl;
        int n(1000);
        //double cnt0(0),cnt1(0),cnt2(0);
        while(keep_going && cnt < 100){
            /*cnt0=0;
            cnt1=0;
            cnt2=0;*/
            int1 = 0;
            int2 = 0;
            if (beta > 0){
                for (int i=0; i<n; i++){
                    eprime = (i-n)*1.0/n;
                    e = 1.0 / eprime + emax + 1;

                    if (e > e1){
                        tmp = a * (exp(e - 2 * log(-eprime))  + exp(e - 2 * log(-eprime) + b*(emax-e))*(b*(emax-e)-1));
                    }
                    if (e < e1 && e > e2){

                        tmp = a * (exp(e - 2 * log(-eprime)) * (1 - 2 * pow(kappa_c/kappa_0,b1)) - exp(e - 2 * log(-eprime) + b*(emax-e))*(b*(emax-e-2*b1/b*log(kappa_c / kappa_0))-1));
                    }
                    if (e < e2){

                        tmp = a * exp(e - 2 * log(-eprime)) * pow(1 - pow(kappa_c/kappa_0,b1),2);
                    }

                    if (beta>=1){
                        tmp *= 1.0 / (1 + exp(beta*e)*pow(mu,-beta));
                        int2 += tmp / (1 + pow(mu,beta)*exp(- beta*e)) * 1.0/n;
                    }
                    else if (beta< 1){
                        tmp *= 1.0 / (1 + exp(beta*e)/mu);
                        int2 += tmp / (1 + mu*exp(- beta*e)) * 1.0/n;
                    }

                    int1 += tmp * 1.0 / n;

                }
            }
            else{
                int1 += std::pow(kappa_0,-gamma)*(kappa_0 * std::pow(kappa_c,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_0*kappa_0*mu)) - kappa_c * std::pow(kappa_0,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_0*kappa_c*mu)));
                int1 += std::pow(kappa_c,-gamma)*(kappa_0 * std::pow(kappa_c,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_0*kappa_c*mu)) - kappa_c * std::pow(kappa_0,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_c*kappa_c*mu)));
                int1 /= 2;
                int2 += std::pow(kappa_0,-1-gamma)*(std::pow(kappa_c,-gamma) * hyp2f1c(gamma,- 1.0 / (kappa_c*kappa_0*mu)) - std::pow(kappa_0,-gamma) * hyp2f1c(gamma,- 1.0 / (kappa_0*kappa_0*mu)));
                int2 += std::pow(kappa_c,-1-gamma)*(std::pow(kappa_c,-gamma) * hyp2f1c(gamma,- 1.0 / (kappa_c*kappa_c*mu)) - std::pow(kappa_0,-gamma) * hyp2f1c(gamma,- 1.0 / (kappa_0*kappa_c*mu)));
                int2 *= (gamma -1) / gamma;
                int2 += mu * std::pow(kappa_0,-gamma)*(std::pow(kappa_0,1-gamma)/(1+kappa_0*kappa_0*mu) - std::pow(kappa_c,1-gamma)/(1+kappa_c*kappa_0*mu));
                int2 += mu * std::pow(kappa_c,-gamma)*(std::pow(kappa_0,1-gamma)/(1+kappa_c*kappa_0*mu) - std::pow(kappa_c,1-gamma)/(1+kappa_c*kappa_c*mu));
                int2 /= 2;

                for (int i=1; i<n; i++){
                    kappa_prime = i * 1.0 / n * (kappa_c - kappa_0) + kappa_0;
                    int1 += std::pow(kappa_prime,-gamma)*(kappa_0 * std::pow(kappa_c,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_0*kappa_prime*mu)) - kappa_c * std::pow(kappa_0,gamma) * hyp2f1d(gamma,- 1.0 / (kappa_c*kappa_prime*mu)));
                    int2 += (gamma -1) / gamma * (std::pow(kappa_prime,-1-gamma)*(std::pow(kappa_c, -gamma) * hyp2f1c(gamma,- 1.0 / (kappa_c*kappa_prime*mu)) - std::pow(kappa_0, -gamma) * hyp2f1c(gamma,- 1.0 / (kappa_0*kappa_prime*mu))));
                    int2 += mu * std::pow(kappa_prime,-gamma)*(std::pow(kappa_0,1-gamma)/(1+kappa_prime*kappa_0*mu) - std::pow(kappa_c,1-gamma)/(1+kappa_c*kappa_prime*mu));
                }

                int1 *= curly_c * curly_c / (gamma - 1) * std::pow(kappa_0 * kappa_c, -gamma) * (kappa_c-kappa_0) * N;
                int1 /= n;

                int2 *= (kappa_c-kappa_0);
                int2 /= n;

                int2 *= curly_c * curly_c * N;
                int2 /= mu * mu;
            }
            //cout << "The contribution of the third term is " << int1 - cnt2 << endl;
            //cout << " --------------- " << endl;
            //cout << "mu = " << mu << " int1 = " << int1  << " int2 = " << int2 << endl;
            //Checks if the goal has been reached, if so, exit the loop.
            if(std::abs(int1 - kavg) < 1e-5){
            keep_going = false;
            }
            else{
            //Newtons method calculates the next step.
            if(beta>=1){mu -= mu / beta * (int1 - kavg)/ int2;}
            else if(beta>0){mu -= mu * (int1 - kavg)/ int2;}
            else{mu -= (int1-kavg) / int2;}
            //If overshoots to negative values enforce trying zero.
            if (mu<=0){mu = 1e-10;}
            ++cnt;
            }
        }
    }
    else{
        //Set initial guess of Newtons method depending on beta.
        if (beta >0){
            if (beta>1){
                mu=beta*sin(M_PI/beta)/(2*M_PI*kavg);
                xi= std::pow(N / (2.0*kavg*kavg),beta); 
            }
            else if (beta == 1){
                mu=1./(2*kavg*log(N));
                xi= N / (2.0*kavg*kavg); 
            }
            else{
                mu=pow(N,beta-1)*(1-beta)/(pow(2,beta)*kavg);
                xi= std::pow(N / 2.0,beta) * 1.0 / (kavg * kavg); 
            }
            while(keep_going && cnt < 100){
                int1 = N * hyp2f1a(beta, - xi * std::pow(mu,-std::max(1.0,beta)));
                int2 = - 1.0 / (mu*std::min(1.0,beta)) * N * (1.0 / (1 + xi * std::pow(mu,-std::max(1.0,beta))) - hyp2f1a(beta,- xi * std::pow(mu,-std::max(1.0,beta))));
                        
                //Checks if the goal has been reached, if so, exit the loop.
                if(std::abs(int1 - kavg) < 1e-5){
                    keep_going = false;
                }
                else{
                    //Newtons method calculates the next step.
                    mu -= (int1 - kavg) / int2;
                    //If overshoots to negative values enforce trying zero.
                    if (mu<=0){mu = 1e-10;}
                    ++cnt;
                }
            }
        }
        else{
            mu = 1.0 / (kavg * (N-kavg));
        }
    }
}

void S1_realization::Fill_thetas(){
    int seed1 = seed + 1;
    mt19937 gen{static_cast<uint32_t>(seed1)};
    uniform_real_distribution<> theta_distrib{0,2*M_PI};
    thetas.clear();
    for (int i=0; i<N; i++){
        thetas.push_back(theta_distrib(gen));
    }
    //sort(thetas.begin(),thetas.end());
}

void S1_realization::Random_Connect(){
    int seed2 = seed + 2;
    mt19937 gen{static_cast<uint32_t>(seed2)};
    uniform_real_distribution<> p_distrib(0,1);

    edges.clear();
    degree.clear();

    edges.resize(2);
    degree.resize(N);

    for (int i=1; i < N; i++){
        /*if (i%10000==1){
            cout << i*1./N << endl;
        }*/

        for (int j=0; j < i; j++){

            if (beta>1){
                xij=N*(M_PI-abs(M_PI-abs(thetas.at(i)-thetas.at(j))))/(2*M_PI*mu*kappas.at(i)*kappas.at(j));
                pij=1/(1+pow(xij,beta));
            }else if (beta>0){
                xij=N*(M_PI-abs(M_PI-abs(thetas.at(i)-thetas.at(j))))/(2*M_PI);
                pij=1/(1+pow(xij,beta)/(mu*kappas.at(i)*kappas.at(j)));
            }else{
                pij=1/(1+1.0/(mu*kappas.at(i)*kappas.at(j)));
            }


            random=p_distrib(gen);
            if (random<pij){
                edges.at(0).push_back(i+1);
                edges.at(1).push_back(j+1);
                degree.at(i)++;
                degree.at(j)++;
            }
        }
    }
}