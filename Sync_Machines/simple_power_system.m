clear;
clc;

% You are given the following power system
% G----|--T--|-----|---->
% A synchronous generator, a step-up transformer, a line, and a load.

% givens
% generator
Sg = 20e6;  % VA
Vg = 6.6e3; % V
Xd = 0.6;   % pu
Xq = 0.5;   % pu
% transformer
St = 30e6;  % VA
Vt1 = 7e3;  % V
Vt2 = 20e3; % V
Vsc = 0.08; % pu
% line
length = 3;                         % km
line_resistance_per_km = 5;
line_reactance_per_km = 12;
line_impedence_per_km = line_resistance_per_km + 1i * line_reactance_per_km; % ohms / km
% load
Sl = 12e6;           % VA   
power_factor = 0.8;  % cos(phi)
lagging = false;     % true for lagging load and false for leading load
Vl = 23e3;           % V

[I_load, I_generator, id, Vt, Ea, Eq, Ei] = simple_sync_machine(Sg, Vg, Xd, Xq, St, Vt1, Vt2, Vsc, length, line_resistance_per_km, line_reactance_per_km, line_impedence_per_km, Sl, power_factor, lagging, Vl);