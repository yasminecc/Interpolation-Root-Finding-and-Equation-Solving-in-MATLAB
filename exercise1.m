% Exercise 1

a = 0;
b = 2;

% Humps function
f = @(x) 1 ./ ((x - 0.3).^2 + 0.01) + 1 ./ ((x - 0.9).^2 + 0.04) - 6;

% True value using integral
true_value = integral(f, a, b);

%% Part (a): Gauss-Legendre Quadrature with n=3

% nodes and weights for [-1, 1]
x = [-sqrt(3/5), 0, sqrt(3/5)];
w = [5/9, 8/9, 5/9];

% nodes for the interval [a, b]
xi = ((b - a) / 2) * x + (a + b) / 2;
wi = ((b - a) / 2) * w;

% Compute the integral estimate
I_gauss = sum(wi .* f(xi));

% true relative error
true_relative_error_gauss = abs((I_gauss - true_value) / true_value) * 100;

fprintf('Integral estimate using Gauss-Legendre Quadrature (n=3): %.6f\n', I_gauss);
fprintf('True relative error Gauss-Legendre: %.6f%%\n', true_relative_error_gauss);

%% Part (b): Dartboard Monte Carlo Integration with n=10000 points

n = 10000;

% Generated random x values within [a, b]
x = a + (b - a) * rand(1, n);

% Evaluate f(x) at these points
f_values = f(x);

% Find approximate min and max values of f(x) over [a, b]
x_grid = linspace(a, b, 10000);
f_grid = f(x_grid);
f_min = min(f_grid);
f_max = max(f_grid);

f_shift = f_values - f_min;
f_shift_max = f_max - f_min;

% Generate random y values between 0 and f_shifted_max
y_mc = f_shift_max * rand(1, n);

% Count the number of points under the shifted curve
num_points = sum(y_mc <= f_shift);

% Estimate the integral of the shifted function
I_shift = (b - a) * f_shift_max * (num_points / n);

% Adjust the integral estimate to get the integral of f(x)
I_mc = I_shift + f_min * (b - a);

% Calculate true relative error for Monte Carlo Integration
true_relative_error_monte = abs((I_mc - true_value) / true_value) * 100;

% Display the result
fprintf('Integral estimate using Dartboard Monte Carlo Integration (n=%d): %.6f\n', n, I_mc);
fprintf('True relative error (Monte Carlo): %.6f%%\n', true_relative_error_monte);