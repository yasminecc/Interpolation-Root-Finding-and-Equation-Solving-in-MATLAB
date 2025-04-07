% Exercise 3

% Parameters
w = 1; d = 1; h = 1;           
xc = 0.5; yc = 0.5; zc = 0.5;  
exact_mass = 0.884323;        

% Density function
rho = @(x, y, z) exp(-((x - xc).^2 + (y - yc).^2 + (z - zc).^2) / 2);

% Grid sizes
grid_sizes = [2, 3, 4];  % For 2x2x2, 3x3x3, and 4x4x4 grids
errors = zeros(size(grid_sizes));
N = grid_sizes.^3;  % Total number of integration points

% calculate the mass for each grid size
for i = 1:length(grid_sizes)
    n = grid_sizes(i);  

    % Gauss-Legendre nodes and weights for each dimension
    [x_nodes, x_weights] = gauss_legendre(n, 0, w);
    [y_nodes, y_weights] = gauss_legendre(n, 0, d);
    [z_nodes, z_weights] = gauss_legendre(n, 0, h);

    % Create a grid of all combinations 
    [X, Y, Z] = ndgrid(x_nodes, y_nodes, z_nodes);
    [Wx, Wy, Wz] = ndgrid(x_weights, y_weights, z_weights);

    % Density at all grid points
    Rho = rho(X, Y, Z);

    % Computing the mass
    mass = sum(Rho(:) .* Wx(:) .* Wy(:) .* Wz(:));

    % Relative true error 
    errors(i) = abs((mass - exact_mass) / exact_mass) * 100;

    % Display the mass for each grid size
    fprintf('Mass estimate for %dx%dx%d grid: %.6f\n', n, n, n, mass);
end

% Plot the relative true error
figure;
plot(N, errors, '-o');
xlabel('Number of integration points (N)');
ylabel('Relative True Error (%)');
title('Relative True Error vs. Number of Integration Points');
grid on;

%% Gauss-Legendre Function

function [points, weights] = gauss_legendre(n, a, b)
  
    switch n
        case 2
            x_standard = [-sqrt(1/3), sqrt(1/3)];
            w_standard = [1, 1];
        case 3
            x_standard = [-sqrt(3/5), 0, sqrt(3/5)];
            w_standard = [5/9, 8/9, 5/9];
        case 4
            x_standard = [-0.861136, -0.339981, 0.339981, 0.861136];
            w_standard = [0.347855, 0.652145, 0.652145, 0.347855];
        otherwise
            error('This function currently supports n=2, 3, or 4 only.');
    end

    % Scale points to the interval [a, b]
    points = ((b - a) / 2) * x_standard + (a + b) / 2;
    weights = ((b - a) / 2) * w_standard;
end
