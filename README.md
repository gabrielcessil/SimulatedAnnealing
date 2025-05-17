# Simulated Annealing Initialization Tool

This repository provides a set of MATLAB functions designed to assist in the initialization of **Simulated Annealing (SA)** algorithms by calculating a suitable starting temperature `T₀`. This is crucial for tuning the algorithm to effectively explore the search space and avoid premature convergence.

The method implemented here is based on an energy-difference-based approach, where the initial temperature `T₀` is estimated using a reference cost difference (`ΔE`) and an initial acceptance probability (`E₀`).

## Purpose

In Simulated Annealing, the choice of the initial temperature `T₀` significantly affects convergence speed and solution quality. The method used here computes `T₀` such that:

```
T₀ = -ΔE / log(E₀)
```

This formula helps define a probabilistic acceptance rule for worse solutions early in the search process, improving the exploration capacity of the SA algorithm.

---

## File Descriptions

- **`calculate_T0.m`**  
  Computes the initial temperature `T₀` based on a user-defined energy difference and initial acceptance probability:
  ```matlab
  T0 = calculate_T0(DELTA_E, E0)
  ```

- **`E0_function.m`**  
  A sample function that evaluates the acceptance probability metric `E₀`. Currently, it returns the maximum value in the input vector:
  ```matlab
  E0 = E0_function(x)
  ```

- **`griewank_function.m`**  
  Implements the **Griewank function**, a common test benchmark for optimization algorithms with a known global minimum at the origin. Useful for evaluating the performance of optimization routines like Simulated Annealing:
  ```matlab
  f = griewank_function(x)
  ```

---

## Example Usage

```matlab
% Example input vector for evaluation
x = rand(1, 10); % Random 10-dimensional vector

% Step 1: Compute E0 using the E0 metric
E0 = E0_function(x);

% Step 2: Define a typical energy difference (e.g., estimated from samples)
DELTA_E = 10;

% Step 3: Calculate initial temperature
T0 = calculate_T0(DELTA_E, E0);

% Step 4: Optionally evaluate the objective function
fitness = griewank_function(x);

% Display results
fprintf('E0 = %.4f, T0 = %.4f, Fitness = %.4f\n', E0, T0, fitness);
```

---

## Requirements

- **MATLAB R2013a** or later (developed and tested on R2013a)

---

## License

This project is provided for educational and research purposes. No specific license has been applied.
