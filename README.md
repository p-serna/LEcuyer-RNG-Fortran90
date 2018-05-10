# LEcuyer-RNG-Fortran90
A version of the Random Number Generator proposed by L'Ecuyer in F90 ready to use. Some comments are in spanish. A quick test:

```
gfortran -c rngstream.f90
gfortran -fopenmp example.f90 rngstream.o -o example
./example
```
