# Harmonic Oscillator on a Lattice
### Project Overview

This project presents a non-perturbative approach to studying a quantum system using numerical Monte Carlo techniques for evaluating the path integral. The focus is on a one-dimensional harmonic oscillator, a system that can be exactly solved analytically, allowing for comparison with our numerical results. The methods used here, however, are general and can be applied to more complex potentials or even quantum field theories with appropriate modifications.
### Features

   - Path Integral Formulation: The quantum system is studied using the path integral formulation, where time is discretized on a lattice, and the integral is computed over finite dimensions using Monte Carlo methods.
   - Monte Carlo Method: The project implements Monte Carlo techniques, specifically the Metropolis algorithm, to evaluate the path integral.
   - Correlation Functions: Calculation of two-point correlation functions and extraction of physical quantities such as energy differences between states and matrix elements.
   - Thermalization and Autocorrelation Time: Analysis of thermalization time and autocorrelation time to ensure accurate sampling.
   - Exact Lattice Solution: Comparison of numerical results with the exact solution on a lattice to validate the methodology.

### Results

The project provides a detailed analysis of the simulation results, including:

   - Selection of parameters and their impact on acceptance rates.
   - Thermalization times as a function of simulation parameters.
   - Calculation of energy levels and matrix elements from the correlation functions.
   - Comparison of numerical results with exact solutions.

### Folders

- <b>devel</b>:            directories used for developing and testing the various 
                 modules. Test programs should be made available here that 
                 can be rerun if so desired

- <b>doc</b>:              collection of notes and other documentation related to the
                 programs in this repository 

- <b>include</b>:          all include files. Typically there is one include file
                 per directory in the modules directory

- <b>main</b>:             collection of main programs

- <b>modules</b>:          source code of all modules 
