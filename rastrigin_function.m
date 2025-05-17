function f = rastrigin_function(x)
    % Rastrigin function
    % Global minimum at x = [0,0,...,0], f = 0
    f = 10 * length(x) + sum(x.^2 - 10 * cos(2 * pi * x));
end