function f = ackley_function(x)
    % Ackley function
    % Global minimum at x = [0,0,...,0], f = 0
    a = 20;
    b = 0.2;
    c = 2 * pi;
    n = length(x);
    f = -a * exp(-b * sqrt(sum(x.^2) / n)) - exp(sum(cos(c * x)) / n) + a + exp(1);
end
