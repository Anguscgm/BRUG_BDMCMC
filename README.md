# BRUG_BDMCMC

To be updated soon! Sorry for the inconvenience!

## Overview
This pure C++ program is written to execute the algorithm in Bayesian Graphical Regression with Birth-Death Markov Process (BGR-BDMCMC) by Yuen et al.
It is constructed and modified above the Gausian Graphical Model Double Metropolis Hastings Birth-Death Markov Process program(ggm_DMH_bdmcmc_ma).
## Installation
BRUG_BDMCMC does not require any installation. It only requires the compiler g++ and the dependecies BLAS and LAPACK for C to be available.
## Compilation and arguments
First, compile the code:
```
g++ BRUG_BDMCMC.cpp -llapack -lblas -std=c++11 -o BRUG_BDMCMC
```
Next, run the compiled code:
```
# For some devices, it may be necessary to update permissions.
chmod 777 ./BRUG_BDMCMC
# Create directory to store outputs.
mkdir demo_out
# Run the program
./BRUG_BDMCMC -i 5000 -b 2500 -p 5 -L 4 -n demo_data/demo_n.txt -B 7 -g demo_data/demo_g_prior.txt -K demo_data/demo_K.txt \
-S demo_data/demo_S.txt -D demo_data/demo_D.txt -z demo_data/demo_beta0.txt -o demo_data/demo_beta1.txt \
-X demo_data/demo_X.txt -e .00000001 -s 123 -O demo_out/
```
The meaning of the arguments are as follows:
```
-i  Integer, the total number of iterations desired.
-b  Integer, the number of burn-in iterations desired.
-p  Integer, the dimension of the covariance matrix.
-L  Integer, the number of groups to be analyzed.
-n  String, the path to the file of the sample sizes by group.
-B  Integer, the initial value of b.
-g  String, the path to the file of g prior.
-K  String, the path to the file of initial K.
-S  String, the path to the file of S.
-D  String, the path to the file of initial D.
-z  String, the path to the file of \beta_0.
-o  String, the path to the file of \beta_1.
-X  String, the path to the file of X, the age of each group.
-e  Double, the value of the threshold.
-s  Integer, the value of the seed.
```
## Demonstration
Please see the ./demo folder for the demo data
## Reference
Mohammadi R, Wit EC (2019). “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software, 89(3), 1–30. doi: 10.18637/jss.v089.i03.
