function [bus_voltages, bus_angles, slack_power] = vectorized_newton(buses, voltages, admittances, assumed_activePower, assumed_reactivePower, max_iterations, tolerance)

arguments (Input)
    buses;
    voltages;
    admittances;
    assumed_activePower;
    assumed_reactivePower;
    max_iterations;
    tolerance;
end

arguments (Output)
    bus_voltages;
    bus_angles;
    slack_power;
end

% Initialize all buses with:
%   voltage = 1 pu
%   angle = 0°
%bus_voltages = ones(length(buses), 1);
bus_angles = zeros(length(buses), 1);

% PV Buses have known voltages
%bus_voltages(buses == 1) = pv_voltages;

bus_voltages = voltages;


i = 0;
while (i < max_iterations)
    P_calc = zeros(length(buses), 1);
    Q_calc = zeros(length(buses), 1);

    for k = 1:length(buses)
        for m = 1:length(buses)
            P_calc(k) = P_calc(k) + bus_voltages(k) * bus_voltages(m) * (real(admittances(k, m)) * cos(bus_angles(k) - bus_angles(m)) + imag(admittances(k, m)) * sin(bus_angles(k) - bus_angles(m)));
            Q_calc(k) = Q_calc(k) + bus_voltages(k) * bus_voltages(m) * (real(admittances(k, m)) * sin(bus_angles(k) - bus_angles(m)) - imag(admittances(k, m)) * cos(bus_angles(k) - bus_angles(m)));
        end
    end

    dP = assumed_activePower(buses ~= 3) - P_calc(buses ~= 3);
    dQ = assumed_reactivePower(buses == 1) - Q_calc(buses == 1);

    F = [dP; dQ];

    if (max(abs(F)) < tolerance)
        break;
    end

    H = zeros(length(buses), length(buses));
    N = zeros(length(buses), length(buses));
    J = zeros(length(buses), length(buses));
    L = zeros(length(buses), length(buses));

    % building the full jacobian sub-matrices
    for k = 1:length(buses)
        for m = 1:length(buses)
            if (k == m)
                H(k, k) = -Q_calc(k) - imag(admittances(k, k)) * bus_voltages(k)^2;
                N(k, k) = P_calc(k)/bus_voltages(k) + real(admittances(k,k))*bus_voltages(k);
                J(k, k) = P_calc(k) - real(admittances(k,k))*bus_voltages(k)^2;
                L(k, k) = Q_calc(k)/bus_voltages(k) - imag(admittances(k,k))*bus_voltages(k);
            else
                H(k, m) = bus_voltages(k) * bus_voltages(m) * (real(admittances(k, m)) * sin(bus_angles(k) - bus_angles(m)) - imag(admittances(k, m)) * cos(bus_angles(k) - bus_angles(m)));
                N(k, m) = bus_voltages(k) * (real(admittances(k, m)) * cos(bus_angles(k) - bus_angles(m)) + imag(admittances(k, m)) * sin(bus_angles(k) - bus_angles(m)));
                J(k, m) = -bus_voltages(k)*bus_voltages(m)*(real(admittances(k, m))*cos(bus_angles(k) - bus_angles(m)) + imag(admittances(k, m))*sin(bus_angles(k) - bus_angles(m)));
                L(k, m) = bus_voltages(k)*(real(admittances(k, m))*sin(bus_angles(k) - bus_angles(m)) - imag(admittances(k, m))*cos(bus_angles(k) - bus_angles(m)));
            end
        end
    end

    % reducing the jacobian
    H = H(buses ~= 3, buses ~= 3);
    N = N(buses ~= 3, buses == 1);
    J = J(buses == 1, buses ~= 3);
    L = L(buses == 1, buses == 1);

    pq_buses = bus_voltages(buses == 1);
    Dv = diag(pq_buses);
    N = N * Dv;
    L = L * Dv;

    jacobian = [H N; 
                J L];

    x = jacobian \ F;
    d_angles = x(1:sum(buses ~= 3));
    rel_dV = x(sum(buses ~= 3)+1:end);
    dV = pq_buses .* rel_dV;

    bus_voltages(buses == 1) = bus_voltages(buses == 1) + dV;
    bus_angles(buses ~= 3) = bus_angles(buses ~= 3) + d_angles;

    i = i + 1; % Increment the iteration counter
end

slack_power = [P_calc; Q_calc];

fprintf('\nBus Voltages:\n');
for k = 1:length(buses)
    fprintf('Bus %d: |V| = %.4f pu, angle = %.4f deg\n', k, bus_voltages(k), rad2deg(bus_angles(k)));
end

end