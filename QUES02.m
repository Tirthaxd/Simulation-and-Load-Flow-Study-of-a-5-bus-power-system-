%{ 
Military Institute of Science and Technology (MIST) 
Dept. of Electrical, Electronic and Communication Engineering 
EECE 306: Power System Laboratory 
Open Ended Lab Project 
 
Group: 03 
 
%} 
clc; clear; close all;  
%% Segment 1: Data Initialization 
bus_data = [1 1.06 0.0 0 0;  % Bus data: [Bus, Voltage (p.u.), Angle, P (MW), Q (MVAR)] 
            2 1.00 0.0 20 20; 
            3 1.00 0.0 -45 -15; 
            4 1.00 0.0 -40 -5; 
            5 1.00 0.0 -60 -10]; 
 
line_data = [1 2 0.02 0.06 0.03;  % Line data: [From Bus, To Bus, R, X, B] 
             1 3 0.08 0.24 0.025; 
             2 3 0.06 0.25 0.02; 
             2 4 0.06 0.18 0.02; 
             2 5 0.04 0.12 0.015; 
             3 4 0.01 0.03 0.01; 
             4 5 0.08 0.24 0.025]; 
 
V_base = [10.6; 10; 10; 10; 10];  % Base voltages in kV 
%% Segment 2: Y-bus Matrix Calculation 
n_lines = size(line_data, 1);  % Number of lines 
n_buses = size(bus_data, 1);   % Number of buses 
Y_bus = zeros(n_buses);        % Initialize Y-bus matrix 
 
for i = 1:n_lines 
    from = line_data(i, 1);  % From bus 
    to = line_data(i, 2);    % To bus 
    Z = line_data(i, 3) + 1i * line_data(i, 4);  % Impedance 
    Y = 1 / Z;  % Admittance 
    B = 1i * line_data(i, 5);  % Shunt admittance 
 
    Y_bus(from, to) = Y_bus(from, to) - Y;  % Off-diagonal elements 
    Y_bus(to, from) = Y_bus(from, to);  % Symmetric update 
    Y_bus(from, from) = Y_bus(from, from) + Y + B;  % Diagonal elements 
    Y_bus(to, to) = Y_bus(to, to) + Y + B;  % Diagonal elements 
end 
 
 
 
disp('Y-bus Matrix:'); 
disp(Y_bus);  % Display Y-bus matrix 
  
%% Segment 3: Power Flow Calculation (Gauss-Seidel) 
base_power = 100;  % Base power in MVA 
V = bus_data(:, 2);  % Initial voltage guess 
P = bus_data(:, 4);  % Real power demand 
Q = bus_data(:, 5);  % Reactive power demand 
tolerance = 1e-5;  % Convergence tolerance 
max_iterations = 100;  % Maximum iterations 
error = 1;  % Initialize error 
iteration = 0;  % Initialize iteration count 
 
while error > tolerance && iteration < max_iterations 
    iteration = iteration + 1;  % Increment iteration 
    V_prev = V;  % Store previous voltages 
 
    for i = 2:n_buses  % Skip slack bus 
        sum_term = 0;  % Initialize sum term 
        for j = 1:n_buses 
            if i ~= j 
                sum_term = sum_term + Y_bus(i, j) * V(j);  % Summation 
            end 
        end 
        V(i) = ((P(i) - 1i * Q(i)) / base_power) / conj(V(i)) - sum_term;  % Voltage update 
        V(i) = V(i) / Y_bus(i, i);  % Normalize 
    end 
 
    error = max(abs(V - V_prev));  % Update error 
end 
 
disp('Converged Voltages (p.u.):'); 
 
 
disp(V);  % Display converged voltages 

%% Segment 4: Voltage Magnitudes and Angles 
voltage_magnitude = abs(V);  % Magnitude of voltages 
voltage_angle = angle(V) * 180 / pi;  % Angle of voltages 
 
fprintf('Bus Voltages (p.u.) and Angles (degrees):\n'); 
for i = 1:n_buses 
    fprintf('Bus %d: |V| = %.4f p.u., Angle = %.2f°\n', i, voltage_magnitude(i), voltage_angle(i));  % Print results 
end 

%% Segment 5: Real Voltages in kV 
V_real = voltage_magnitude .* V_base;  % Calculate real voltages in kV 
 
fprintf('\nReal Voltages in kV:\n'); 
for i = 1:n_buses 
    fprintf('Bus %d: Voltage = %.4f kV, Angle = %.2f°\n', i, V_real(i),voltage_angle(i));  % Print results 
end

%% Segment 6: Generation and Load Power at Each Bus 
fprintf('\nGeneration and Load Power at Each Bus:\n'); 
for i = 1:n_buses 
    S_bus = V(i) * conj(Y_bus(i, :) * V);  % Calculate complex power 
    P_bus = real(S_bus) * base_power;  % Real power 
    Q_bus = imag(S_bus) * base_power;  % Reactive power 
     
    if P_bus >= 0 
        fprintf('Bus %d: Generation Power = %.2f MW, Reactive Power = %.2f MVAR\n', i, P_bus, Q_bus); 
    else 
        fprintf('Bus %d: Load Power = %.2f MW, Reactive Power = %.2f MVAR\n', i, P_bus, -Q_bus); 
    end 
end 


