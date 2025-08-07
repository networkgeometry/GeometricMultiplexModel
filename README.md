# Geometric Multiplex Model

## Introduction

This repository can be used to generate an instance of the _Geometric Multiplex  (GMM). This is an extention of the geometric $\mathbb{S}^1$-model with $L$ layers where the correlations at the level of the hidden coordinates between the layers can be tuned. The full model can be found in [Kleineberg2016](https://doi.org/10.1038/nphys3812). 

## Installation

Requirements 

* A C++11 (or newer) compliant compiler
* Boost C++ libraries
    * boost/math/special_functions/lambert_w.hpp
    * boost/math/special_functions/erf.hpp

```
# Unix (Linux / MAC OS)
g++ -O3 -std=c++11 include_else/hyp2f1.cpp include_me/S1_realization.cpp main.cpp -o GENMULTIPLEX
```

## Usage

### Format of input files

There are two input files containing S1 parameters and interlayer correlations respectively.

The file parameters.txt contains $L$ lines, one for each layer. It must have the following structure:

```
[gamma layer 1] [beta layer 1] [average degree layer 1] [random seed layer 1]
[gamma layer 2] [beta layer 2] [average degree layer 2] [random seed layer 2]
...
[gamma layer L] [beta layer L] [average degree layer L] [random seed layer L]
```

***Note:*** This version of the GMM supports two degree distributions:<br>
* Pareto distribution with cut-off $$\rho(\kappa)=\frac{(\gamma-1)\kappa_0^{\gamma-1}}{1-\left(\frac{\kappa_c}{\kappa_0}\right)^{1-\gamma}}\kappa^{-\gamma}.$$ This  distribution is chosen automatically if $\gamma>2$. 
* Uniform distribution $\rho(\kappa) = \delta\left(\kappa-\langle k\rangle\right)$. This distribution is chosen automatically if $\gamma\leq 2$.

The file correlations.txt contains $L-1$ lines, one for each consecutive pair of layers. It must have the following structure:

```
[nu layers 1,2] [g layers 1,2]
[nu layers 2,3] [g layers 2,3]
...
[nu layers L-1,L] [g layers L-1,L]
```

### Running the code

Running `GENMULTIPLEX` can be done with or without flags

```
# Command line with flags
./GENMULTIPLEX -N <network size> -pf <parameter file name> -cf <correlation file name>

# Command line without flags
./GENMULTIPLEX <network size> <parameter file name> <correlation file name>
```

### Output files

The program outputs the multiplex using a set of files. For each layer `i`, a hidden coordinate file `layer<i>.coord` is generated which contains two columns. The first contains the angular coordinates and the second the hidden degrees. Additionally, for each layer `i` an edge file `layer<i>.edge` is created.

## Publications 

Please cite:

_Hidden geometric correlations in real multiplex networks_<br>
Kaj Kolja Kleineberg, Marián Boguñá, M. Ángeles Serrano and Fragkiskos Papadopoulos
<br>
Nature Physics __12__, 1076-1081 (2016)<br>
[Full text](https://doi.org/10.1038/nphys3812) | [arXiv](https://doi.org/10.48550/arXiv.1601.04071)

_Multiplexity amplifies geometry in networks_<br>
Jasper van der Kolk, Dmitri Krioukov, Marián Boguñá and M. Ángeles Serrano<br>
arXiv preprint arXiv:2505.17688<br>
[arXiv](https://doi.org/10.48550/arXiv.2505.17688)
