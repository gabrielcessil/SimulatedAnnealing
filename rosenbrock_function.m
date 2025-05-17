function f = rosenbrock_function(x)
    % Rosenbrock function
    % Global minimum at x = [1,1,...,1], f = 0
    f = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (x(1:end-1) - 1).^2);
end