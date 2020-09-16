c   This code computes the sensitivity of a model to variations (or equivalently uncertainties) of its own parameters. The estimators of variance-based sensitivity analysis are used to achieve this goal. These are the so-called Si, Sti, Sij, Stij where i and j are the ith and jth parameters of the model. In addition, one can compute the mean an standard deviation of the model's output with respect to the parameters variations.
c   Si is the 1st-order sensitivity index of the model with respect to its ith parameter.
c   Sti is the total-order sensitivity index of the model with respect to its ith parameter.
c   Sij is the 2nd-order sensitivity index of the model with respect to its ith and jth parameters.

c   For more information about the definitions used in this code, refer to the Wikipedia page on variance-based sensitivity analysis: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis

c   The sensitivity indices are evaluated by (pseudo-)randomly sampling the parameter phase space. However, rather than using a basic Monte Carlo method, we use low-discrepancy Sobol sequences to simulate the sampling. This allows to cover the phase space in a more uniform fashion than random sampling with a given number of samples, improving the confidence one can have in the estimators.

c   This code provides a main file where a toy model with 3 varying parameters a,b,c and 3 fixed inputs x,y,z. One can easily test the code by modifying the relations between the parameters in the model definition. Adding/removing varying parameters will require modifications of the main file and the 'par_estimators' module.

c   Sensitivity analysis can quickly become computationally expensive when the model is complex and requires long computation times, since the model is run many times. The loop where the model is called can be parallelized using the OpenMP library in order to accelerate the computations.

c   To run the test code with gfortran, use the following commands: gfortran -c par_estimators.for && gfortran main.for par_estimators.o && ./a.out
c   To use OpenMP for code parallelization in the par_estimators module with gfortran, use the following commands: gfortran -fopenmp -c par_estimators.for && gfortran -fopenmp main.for par_estimators.o && ./a.out
c   In order for parallelization to work, one needs to uncomment the '!$omp parallel do' loop instructions in par_estimators, and eventually the 'use omp_lib' instruction in order to use specific OpenMP functionalities.


      program main
      use par_estimators
      implicit none

      integer :: seed
      integer N,d,i
      double precision x,y,z
      double precision meanAll,varAll,acc,resmean,resvar
      double precision V1,S1,V2,S2,V3,S3,Vt1,St1,Vt2,St2,Vt3,St3,V12,S12
     &,V13,S13,V23,S23

c   Generation of normalised random numbers (i.e. between 0 and 1) using a randomized seed (fixed seed => fixed "random" numbers since they come from a pseudorandom sequence) thanks to the CPU clock. Randomized seed allows to re-run the program with different values each time
c      call system_clock(seed)
c      seed=int(seed)
c      write(*,*) 'seed=',seed
c      call srand(seed)

c   Fixed seed to compare results between different runs
      seed=8600586
      write(*,*) 'seed=',seed

c   Inputs x,y,z of the toy model
      x=1d0
      y=1d0
      z=1d0

c   Minimal and maximal values for our uncertain model parameters : a,b,c => d=3
      amin=0.0d0
      amax=1.0d0
      bmin=0.0d0
      bmax=1.0d0
      cmin=0.0d0
      cmax=1.0d0
c   For gaussian-based distribution of parameters, one must provide the mean and variance of the parameter's distribution
c      cm=0.5d0
c      cv=0.25d0

c   Initial number of sample points taken in the normalised phase-space for the analysis, which is also the number of sample points for each increment in the convergence computation of the estimators
      N=10000
c   Dimensions of the phase space (i.e. number of varying model parameters) (for d!=3, other changes need to be made manually in the code as the number of arrays in the Sobol sequence is related to d)
      d=3
c   Accuracy goal for the convergence test (here 1%)
      acc=0.01d0
      write (*,*) 'N=',N,',   d=',d,char(10)

c   Call of the various estimators for a model with d=3 varying parameters, the estimators are defined in the 'par_estimators' module
c   Here we use a toy model for demonstration purposes

c   convMeanVar computes the mean and variance of the model's output for the provided parameter distributions and inputs
      call convMeanVar(N,d,acc,toymodel,resmean,resvar)
      write(*,*) char(10)
c   convEstim computes the desired (non-averaged) estimator of variance-based sensitivity analysis depending on which flag is provided: Vi or Vti for first-order or total-order sensitivity index respectively, as well as mean and variance
      call convEstim(1,N,d,acc,toymodel,'Vi ',V1,S1,resmean,resvar)
      call convEstim(1,N,d,acc,toymodel,'Vti ',Vt1,St1,resmean,resvar)
      call convEstim(2,N,d,acc,toymodel,'Vi ',V2,S2,resmean,resvar)
      call convEstim(2,N,d,acc,toymodel,'Vti ',Vt2,St2,resmean,resvar)
      call convEstim(3,N,d,acc,toymodel,'Vi ',V3,S3,resmean,resvar)
      call convEstim(3,N,d,acc,toymodel,'Vti ',Vt3,St3,resmean,resvar)
c   2nd-order sensitivity indices are computed by combining 1st-order and total-order sensitivity indices. This requires to assume 3rd- and higher-order sensitivity indices are negligible
      S12=Sij(S1,S2,S3,St1,St2,St3)
      S13=Sij(S1,S3,S2,St1,St3,St2)
      S23=Sij(S2,S3,S1,St2,St3,St1)
      write(6,*) '  S12         S13         S23'
      write(6,16) S12,S13,S23

 13   format(' ',9f8.3)
 14   format(' ',f15.3)
 15   format(8(1pe13.5,2x))
 16   format(8(1pe10.2,2x))

      contains

c   Toy model of dimension d=3
      double precision function toymodel()
      toymodel=ap*x+bp*y+cp*z
      end function toymodel

      end program main
