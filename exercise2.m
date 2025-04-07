% Exercise 2

air_density = 1.2;  

% Radius (cm) and Velocity (m/s)
r_cm = [0, 1.6, 3.2, 4.8, 6.4, 7.47, 7.87, 7.95, 8];  
v = [10, 9.69, 9.30, 8.77, 7.95, 6.79, 5.57, 4.89, 0];  

% Convert radius to meters
r = r_cm / 100;

% f(r) = density * v(r) * r
f = air_density * v .* r;

n = length(r);

% Initialize total integral
I_total = 0;

simpson_indices = [1, 3, 5];

for i = simpson_indices
    h = r(i+2) - r(i);
    I_segment = (h / 6) * (f(i) + 4*f(i+1) + f(i+2));
    I_total = I_total + I_segment;
end

%Trapezoidal rule for the last interval [r(8), r(9)]
h_trapz = r(n) - r(n-1);
I_trapz = (h_trapz / 2) * (f(n-1) + f(n));
I_total = I_total + I_trapz;

mass_flow_rate = 2 * pi * I_total;

fprintf('Mass flow rate: %.6f kg/s\n', mass_flow_rate);
