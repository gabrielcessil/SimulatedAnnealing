close all;

% Initialize the SA algorithm
sa = Continuos_SA(  'local_suc', 5, ...
            'local_attempts', 30, ...
            'max_it', 10000, ...
            'temp_dec', 0.2, ...
            'Wanted_Cost', -1, ...
            'Max_Cost', 1e6, ...
            'verbose', false, ...
            'convergence_threshold', 1e-6, ...
            'random_seed', 2021);

% Run the optimization for sphere function

initial_solution = [5, -5];


fprintf('Initial Solution for sphere_function: %s with cost %.6f\n', mat2str(initial_solution), sphere_function(initial_solution));
best_solution = sa.run(initial_solution, @sphere_function, 1);
fprintf('Best Solution found for sphere_function: %s with cost %.6f\n\n', mat2str(best_solution), sphere_function(best_solution));


fprintf('Initial Solution for ackley_function: %s with cost %.6f\n', mat2str(initial_solution), ackley_function(initial_solution));
best_solution = sa.run(initial_solution, @ackley_function, 1);
fprintf('Best Solution found for ackley_function: %s with cost %.6f\n\n', mat2str(best_solution), ackley_function(best_solution));


fprintf('Initial Solution for griewank_function: %s with cost %.6f\n', mat2str(initial_solution), griewank_function(initial_solution));
best_solution = sa.run(initial_solution, @griewank_function, 1);
fprintf('Best Solution found for griewank_function: %s with cost %.6f\n\n', mat2str(best_solution), griewank_function(best_solution));


fprintf('Initial Solution for rastrigin_function: %s with cost %.6f\n', mat2str(initial_solution), rastrigin_function(initial_solution));
best_solution = sa.run(initial_solution, @rastrigin_function, 1);
fprintf('Best Solution found for rastrigin_function: %s with cost %.6f\n\n', mat2str(best_solution), rastrigin_function(best_solution));

fprintf('Initial Solution for rosenbrock_function: %s with cost %.6f\n', mat2str(initial_solution), rosenbrock_function(initial_solution));
best_solution = sa.run(initial_solution, @rosenbrock_function, 1);
fprintf('Best Solution found for rosenbrock_function: %s with cost %.6f\n\n', mat2str(best_solution), rosenbrock_function(best_solution));
