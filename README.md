# Var_sensi_analysis: A FORTRAN code for variance-based sensitivity analysis

## Goal of the project

This code computes the sensitivity of a model to variations (or equivalently uncertainties) of its own parameters. The estimators of variance-based sensitivity analysis are used to achieve this goal. These are the so-called *Si*, *Sti*, *Sij* where *i* and *j* are the *i*th and *j*th parameters of the model. In addition, one can compute the mean an standard deviation of the model's output with respect to the parameters variations.
* *Si* is the 1st-order sensitivity index of the model with respect to its *i*th parameter.
* *Sti* is the total-order sensitivity index of the model with respect to its *i*th parameter.
* *Sij* is the 2nd-order sensitivity index of the model with respect to its *i*th and *j*th parameters.

For more information about the definitions used in this code, refer to [the Wikipedia page on variance-based sensitivity analysis](https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis).

The sensitivity indices are evaluated by (pseudo-)randomly sampling the parameter phase space. However, rather than using a basic Monte Carlo method, we use low-discrepancy Sobol sequences to simulate the sampling. This allows to cover the phase space in a more uniform fashion than random sampling with a given number of samples, improving the confidence one can have in the estimators. Such approaches are part of quasi-Monte Carlo methods.

## Contents

The project contains two files:
* [main.for](./main.for) is the main file where a toy model is defined to make tests with and that calls the subroutines computing the estimators
* [par_estimators.for](./par_estimators.for) is a module where the subroutines that compute estimators and the related helper functions are defined

## Getting started

This code provides a main file where a toy model with 3 varying parameters a,b,c and 3 fixed inputs x,y,z. One can easily test the code by modifying the relations between the parameters in the model definition. Adding/removing varying parameters will require modifications of both the main file and the 'par_estimators' module.

To run the test code with gfortran, use the following commands: `gfortran -c par_estimators.for && gfortran main.for par_estimators.o && ./a.out`

### Parallel computing with OpenMP

Sensitivity analysis can quickly become computationally expensive when the model is complex and requires long computation times, since the model is run many times. The loop where the model is called can be parallelized using the OpenMP library in order to accelerate the computations.

In order for parallelization to work, one needs to uncomment the `!$omp parallel do` loop instructions in [par_estimators.for](./par_estimators.for), and eventually the `use omp_lib` instruction in order to use specific OpenMP functionalities.

To use OpenMP for code parallelization in the par_estimators module with gfortran, use the following commands: `gfortran -fopenmp -c par_estimators.for && gfortran -fopenmp main.for par_estimators.o && ./a.out`
