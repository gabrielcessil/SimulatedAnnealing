classdef Continuos_SA
    % Continuos_SA: Continuous Simulated Annealing optimization algorithm.
    % This class implements a simulated annealing algorithm for minimizing
    % a cost function in a continuous search space.

    properties
        local_suc % Number of local
        local_attempts % Number of local attempts per global iteration
        max_it % Maximum number of global iterations
        alpha % Cooling rate (temperature decay factor)
        T0 % Initial temperature
        Wanted_Cost % Target cost (stopping criterion)
        Max_Cost % Maximum acceptable cost
        verbose % Display flag for debugging information
        convergence_threshold % Threshold for cost change to stop optimization
        random_seed % Random seed for reproducibility
        E0 % Customizable function for calculating initial temperature
        scale % Parameter's perturbation scale
        FO_h % History of objective function values
        T_h % History of temperatures
        sol_h % History of solutions
        Perturbationcao_h % History of Perturbationtions
        avalia % Function handle for cost evaluation
    end

    methods
        function obj = Continuos_SA(varargin)
            % Constructor for Continuos_SA class.
            % Inputs:
            %   varargin: Name-value pairs for parameters

            % Create input parser
            p = inputParser;

            % Define name and default value for each parameter
            %(REPLACE WITH addParameter IN NEWER MATLAB VERSIONS)
            addParamValue(p, 'local_suc', 3);
            addParamValue(p, 'local_attempts', 10);
            addParamValue(p, 'max_it', 100);
            addParamValue(p, 'temp_dec', 0.2);
            addParamValue(p, 't0', -1);
            addParamValue(p, 'Wanted_Cost', -1);
            addParamValue(p, 'Max_Cost', sqrt(realmax) - 1);
            addParamValue(p, 'verbose', false);
            addParamValue(p, 'convergence_threshold', 1e-6);
            addParamValue(p, 'random_seed', 2021);
            addParamValue(p, 'E0', 0.8);

            % Parse inputs
            parse(p, varargin{:});
            params = p.Results;

            % Assign properties
            obj.local_suc = params.local_suc;
            obj.local_attempts = params.local_attempts;
            obj.max_it = params.max_it;
            obj.alpha = params.temp_dec;
            obj.T0 = params.t0;
            obj.Wanted_Cost = params.Wanted_Cost;
            obj.Max_Cost = params.Max_Cost;
            obj.verbose = params.verbose;
            obj.convergence_threshold = params.convergence_threshold;
            obj.random_seed = params.random_seed;
            obj.E0 = params.E0;
            
            obj.scale = NaN;
            % Set random seed
            rng(obj.random_seed);
            
        end
        
        function DELTA_E = Delta_E(obj, Xinit, Fcusto, Imax)
            positive_counter = 0; % Count positive cost increments
            DELTAS = zeros(1,Imax); % Initialize saving array
            
            X_old = Xinit; 
            custo_old = Fcusto(X_old);
            for i = 1:Imax
                
                if obj.verbose
                    fprintf(' Calculating DELTA_E: %d/%d\n', i, Imax);
                end
                 
                X_new = obj.Perturbation(X_old);
                
                try
                    custo_new = Fcusto(X_new);
                catch
                    % Ignore new solution in case of error                      
                    % (prohibitive new solution)
                    custo_new = custo_old;
                    continue
                end

                % If the cost increased, count transition 
                if(custo_new - custo_old > 0 )
                    DELTAS(i) = custo_new - custo_old;
                    positive_counter = positive_counter+1;
                end
                
                % Update info for next iteration
                X_old = X_new; 
                custo_old = custo_new;

            end

            % Average across the random perturbations with cost increase
            DELTA_E = sum(DELTAS)/positive_counter;  
        end
        
        
        function T0 = Temp_inic(obj, E0, s0, Cost_evaluate, I)
            delta_E = obj.Delta_E(s0, Cost_evaluate, I);
            T0 = - delta_E / log(E0);
        end

        function [aux, dif] = Perturbation(obj, s)
            aux = s; % No deepcopy needed for basic types
            i = randi(length(s)); % Randomly select an index to perturb
            aux(i) = s(i) + randn() * abs(obj.scale); % Normal distribution
            dif = abs(s(i) - aux(i));
        end
        
        function plotHistory(obj)
            figure;
            subplot(4, 1, 1); plot(obj.FO_h); title('Objective Function (FO)');
            subplot(4, 1, 2); plot(obj.sol_h); title('Solution History');
            subplot(4, 1, 3); plot(obj.T_h); title('Temperature History');
            subplot(4, 1, 4); plot(obj.Perturbationcao_h); title('Perturbationtion History');
            
            % After creating your figure, add:
            print('history', '-dpng', '-r300'); % Saves as 300 DPI PNG0
        end
        function plotSolutionScatter(obj, Cost_evaluate, color_style)
            % Creates visualization with selectable coloring style
            % color_style: 'background' (default) or 'markers'

            % Set default style if not provided
            if nargin < 3
                color_style = 'background';
            end

            % Check if there's any data to plot
            if isempty(obj.sol_h)
                warning('No solution history data available to plot.');
                return;
            end

            % Extract solution data
            solutions = obj.sol_h;
            costs = obj.FO_h;
            num_vars = size(solutions, 2);

            % Create figure
            figure;

            % For each pair of consecutive variables
            for i = 1:min(num_vars-1, 3) % Show max 3 parameter pairs
                subplot(1, min(num_vars-1, 3), i);
                hold on;

                % Get current variable pair range with padding
                x_var = solutions(:, i);
                y_var = solutions(:, i+1);
                padding = 0.1;
                x_pad = padding*(max(x_var)-min(x_var));
                y_pad = padding*(max(y_var)-min(y_var));
                x_vals = linspace(min(x_var)-x_pad, max(x_var)+x_pad, 50);
                y_vals = linspace(min(y_var)-y_pad, max(y_var)+y_pad, 50);

                if strcmpi(color_style, 'background')
                    % BACKGROUND COLORING STYLE
                    [X, Y] = meshgrid(x_vals, y_vals);

                    % Calculate costs for grid
                    Z = zeros(size(X));
                    base_point = mean(solutions, 1);
                    for row = 1:size(X,1)
                        for col = 1:size(X,2)
                            test_point = base_point;
                            test_point(i) = X(row,col);
                            test_point(i+1) = Y(row,col);
                            Z(row,col) = Cost_evaluate(test_point);
                        end
                    end

                    % Plot cost surface
                    contourf(X, Y, Z, 20, 'LineColor', 'none');
                    caxis([min(costs) max(costs)]); % Consistent color scaling

                    % Plot markers (white with black borders)
                    marker_face = 'w';
                else
                    % MARKER-ONLY COLORING STYLE
                    % Just plot the axes
                    xlim([min(x_var)-x_pad max(x_var)+x_pad]);
                    ylim([min(y_var)-y_pad max(y_var)+y_pad]);

                    % Plot markers (colored by cost)
                    marker_face = 'flat';
                end

                % Plot all solution points
                scatter(solutions(:,i), solutions(:,i+1), 60, costs, ...
                        'MarkerFaceColor', marker_face, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 1);

                % Highlight special points
                scatter(solutions(1,i), solutions(1,i+1), 150, ...
                        'MarkerFaceColor', marker_face, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 2, ...
                        'Marker', 'square');

                scatter(solutions(end,i), solutions(end,i+1), 150, ...
                        'MarkerFaceColor', marker_face, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 2, ...
                        'Marker', 'd');

                % Labels and colorbar
                xlabel(['Parameter ', num2str(i)]);
                ylabel(['Parameter ', num2str(i+1)]);
                title(['Params ', num2str(i), '-', num2str(i+1)]);
                grid on;

                if strcmpi(color_style, 'background')
                    colorbar;
                end
            end

            % Set colormap and figure size
            colormap(jet);
            set(gcf, 'Position', [100 100 800*min(num_vars-1, 3) 600]);

            % Add legend
            if strcmpi(color_style, 'background')
                legend('Cost landscape', 'Solutions', 'Start', 'End');
            else
                legend('Solutions', 'Start', 'End');
            end
            
            print('Map', '-dpng', '-r300'); % Saves as 300 DPI PNG0
        end
        function best_solution = run(obj, s0, Cost_evaluate, perturbation_scale)
            % Run the simulated annealing algorithm
            
            % Initialize problem scale
            obj.scale = perturbation_scale;
            
            % Initialize history arrays with preallocation
            obj.FO_h = zeros(1, obj.max_it * obj.local_attempts); % Preallocate for worst-case scenario
            obj.T_h = zeros(1, obj.max_it);
            obj.sol_h = zeros(obj.max_it * obj.local_attempts * length(obj.local_suc), length(s0)); % Preallocate for worst-case scenario
            obj.Perturbationcao_h = zeros(1, obj.max_it * obj.local_attempts); % Preallocate for worst-case scenario
            
            
            
            % Initialize temperature
            if obj.T0 <= 0
                T = obj.Temp_inic(obj.E0, s0, Cost_evaluate, 30);
            else
                T = obj.T0;
            end

            Solution = s0;
            current_cost = Cost_evaluate(Solution);
            best_cost = current_cost;
            best_solution = Solution;
            history_index = 1; % Index for preallocated history arrays
            
            % Operate M different locals
            for j = 1:obj.max_it % Iterate up to the defined max number of iterations
                nsucc = 0; % Reset local success counter
                if obj.verbose
                    fprintf('New Location: %i /%i \n',j,obj.max_it);
                end
                
                % Local Search
                for i = 1:obj.local_attempts
                    [new_Solution, dif] = obj.Perturbation(Solution);
                 
                    new_cost = Cost_evaluate(new_Solution);
                 
                    cost_increment = new_cost - current_cost;

                    if new_cost > obj.Max_Cost
                        if obj.verbose
                            fprintf('Cost above the tolerance: %f\n', obj.Max_Cost);
                        end
                        continue
                    end
                    
                    
                    % ACCEPTANCE FUNCTION
                    acceptance_probability = min(1, exp(-cost_increment / T)); 
                    random_value = rand();
                    
                    
                    if obj.verbose
                        fprintf('Solution %s perturbed to %s.\n', mat2str(Solution), mat2str(new_Solution));
                        fprintf('-- Solution parameter step= %.2f\n', dif);
                        fprintf('-- Cost= %.2f\n', new_cost);
                        fprintf('-- Cost Decay= %.2f\n', cost_increment);
                        fprintf('-- T= %.2f\n', T);
                        fprintf('-- Acceptance probability= %.6f\n', acceptance_probability);
                        fprintf('-- Random value= %.8f\n', random_value);
                    end

                    % Success case: cost decay or stochastic acceptance
                    if acceptance_probability >= random_value
                        fprintf('Success at iteration %d: cost=%.4f\n', history_index, new_cost);
                        % Update solution informations
                        obj.FO_h(history_index) = new_cost;
                        obj.sol_h(history_index, :) = new_Solution;
                        obj.Perturbationcao_h(history_index) = dif;    
                        obj.T_h(history_index) = T;
                        
                        % COOLING SCHEDULE
                        T = max(best_cost/log(j+1), 0.01);
                        fprintf('T %.2f = best_cost %.2f / log(%d) %.2f\n',T,best_cost,j, log(j));
                        
                        if new_cost < best_cost
                            best_cost = new_cost;
                            best_solution = new_Solution;
                        end
                        
                        Solution = new_Solution;
                        current_cost = new_cost;
                        nsucc = nsucc + 1;
                        
                        history_index = history_index + 1;
                
                    end
                    
                    % Break local search
                    if nsucc >= obj.local_suc
                        if obj.verbose
                            fprintf('Enough local sucesses. Going to next location.');
                        end
                        break; % Stop if enough local successes
                    end
                end
               
                
                % Check global convergence
                if j > 2 && abs(obj.FO_h(history_index - 1) - obj.FO_h(history_index - 2)) < obj.convergence_threshold
                    if obj.verbose
                        fprintf('Convergence reached at iteration %d\n. Stopping search.', j);
                    end
                    break;
                end
               
                
                
                if nsucc == 0
                    if obj.verbose
                        fprintf('No solution suceeded on iteration %d. Stopping search\n', j);
                    end
                    break;
                end

                
                if current_cost <= obj.Wanted_Cost
                    if obj.verbose
                        fprintf('Desired cost achieved on iteration %d\n. Stopping search', j);
                    end
                    break;
                end
            end
            % Trim unused preallocated space (important!)
            obj.FO_h = obj.FO_h(1:history_index - 1);
            obj.T_h = obj.T_h(1:history_index - 1);
            obj.sol_h = obj.sol_h(1:history_index - 1, :);
            obj.Perturbationcao_h = obj.Perturbationcao_h(1:history_index - 1);
            
            
            
            % Plots
            obj.plotHistory();
            obj.plotSolutionScatter(Cost_evaluate);




        end

        
    end
end


