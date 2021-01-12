# BRUG_BDMCMC

To be updated and references will be added soon! Sorry for the inconvenience!

## Overview

## Installation
BRUG_BDMCMC does not require any installation. It only requires the compiler g++ and the dependecy LAPACK for C to be available.
## Demonstration
First, compile the code:
```
g++ BRUG_BDMCMC.cpp -llapack -lblas -std=c++11 -o BRUG_BDMCMC
```
Next, run the compiled code:
```
# For some devices, it may be necessary to updatee permissions.
chmod 777 ./BRUG_BDMCMC
# Run the program
./BRUG_BDMCMC -i 5000 -b 2500 -p 4 -L 3 -n demo_n.txt -B 7 -g demo_g_prior.txt -t demo_Ts.txt -T demo_Ti.txt -K demo_K.txt \
-d demo_Ds.txt -D demo_D.txt -z demo_beta0.txt -o demo_beta1.txt-X -e .00000001 -s 123
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
-t  String, the path to the file of T_s.
-T  String, the path to the file of T_i.
-K  String, the path to the file of initial K.
-d  String, the path to the file of initial D_s.
-D  String, the path to the file of initial D.
-z  String, the path to the file of \beta_0.
-o  String, the path to the file of \beta_1.
-X  String, the path to the file of X, the age of each group.
-e  Double, the value of the threshold.
-s  Integer, the value of the seed.
```
## Reference
Mohammadi R, Wit EC (2019). “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software, 89(3), 1–30. doi: 10.18637/jss.v089.i03.
