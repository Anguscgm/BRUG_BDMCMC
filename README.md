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
-i
-b
-p
-L
-n
-B
-g
-t
-T
-K
-d
-D
-z
-o
-X
-e
-s
```
## Reference
Mohammadi R, Wit EC (2019). “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software, 89(3), 1–30. doi: 10.18637/jss.v089.i03.
