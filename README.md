# LEcuyer-RNG-Fortran90

## About

A version in Fortran 90 of a pseudo-random number generator, more precisely, a combined multiple recursive generator (CMRG), MRG32k3a, proposed by [L'Ecuyer](http://www-labs.iro.umontreal.ca/~lecuyer/) in

> L'ecuyer P. Good parameters and implementations for combined multiple recursive random number generators. Operations Research. 1999 Feb;47(1):159-64.

This algorithm is well suited for parallel computing as it can initialize many long streams and substreams, with a total period length of 2<sup>191</sup> and it passes [diehard tests](https://en.wikipedia.org/wiki/Diehard_tests).

## Getting started

File rngstream.f90 is the file for the definition of the module. A quick example can be found in example.f90

```
gfortran -c rngstream.f90
gfortran -fopenmp example.f90 rngstream.o -o example
./example
```
## Disclaimer

Even if it passes diehard tests, it does not mean it suits for any other particular requirements. It seems to work well with Monte Carlo simulations. Some comments are still in spanish.
