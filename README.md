# Simulated Annealing Initialization Tool

This repository provides a set of MATLAB functions designed to assist in the initialization of **Simulated Annealing (SA)** algorithms by calculating a suitable starting temperature `T₀`. This is crucial for tuning the algorithm to effectively explore the search space and avoid premature convergence.

The method implemented here is based on an energy-difference-based approach, where the initial temperature `T₀` is estimated using a reference cost difference (`ΔE`) and an initial acceptance probability (`E₀`).

## Purpose

In Simulated Annealing, the choice of the initial temperature `T₀` significantly affects convergence speed and solution quality. The method used here computes `T₀` such that:

T₀ = -ΔE / log(E₀)


This formula helps define a probabilistic acceptance rule for worse solutions early in the search process, improving the exploration capacity of the SA algorithm.

---

## File Descriptions

- **`calculate_T0.m`**  
  Computes the initial temperature `T₀` based on a user-defined energy difference and initial acceptance probability:
  ```matlab
  T0 = calculate_T0(DELTA_E, E0)
