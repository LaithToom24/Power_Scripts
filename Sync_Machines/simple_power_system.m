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

% unknowns 
% loading_angle
% load_angle
% I      load current
% id     direct axis current
% Vt     transformer voltage
% Ea     terminal voltage
% Eq     quadature voltage
% Ei     internal voltage
% I_pu   
% id_pu 
% Vt_pu  
% Ea_pu  
% Eq_pu  
% Ei_pu  

% preliminary calculations
S_base = St; % set base apparent power to transformer apparent power
             % to avoid recalculation of transformer values
Zb1 = Vt1^2 / S_base;
Zb2 = Vt2^2 / S_base;
line_impedence = length * line_impedence_per_km;
line_impedence_pu = line_impedence / Zb2;
Xt_pu = Vsc;
Sg_pu = Sg / S_base;
Sl_pu = Sl / S_base;
St_pu = St / S_base;
Vg_pu = Vg / Vt1;
Vl_pu = Vl / Vt2;
% recalculating generator reactacnes xq and xd
Zb_g = Vg^2 / Sg;
Xq_new = Xq * (Zb_g / Zb1);
Xd_new = Xd * (Zb_g / Zb1);
I_pu = (Sl_pu / Vl_pu);
load_angle = acos(power_factor);
if (lagging == true)
    load_angle = -load_angle;
end
I_pu = I_pu * exp(1i * load_angle);
Vt_pu = I_pu * line_impedence_pu + Vl_pu;
Ea_pu = I_pu * 1i * Xt_pu + Vt_pu;  
Eq_pu = I_pu * 1i * Xq_new + Ea_pu;
loading_angle = angle(Eq_pu) - angle(Ea_pu);
id_pu = abs(I_pu) * sin(loading_angle + angle(Ea_pu) - load_angle);
Ei_pu = abs(Eq_pu) + (Xd_new - Xq_new) * id_pu;

% converting to SI
Ei = Ei_pu * Vt1; 
Eq = Eq_pu * Vt1;
Ea = Ea_pu * Vt1;
Vt = Vt_pu * Vt2;
I1 = I_pu * S_base/Vt1;
I2 = I_pu * S_base/Vt2;

% confirming power equilibrium
% Check power equilibrium
P_gen_pu = ( abs(Ei_pu)*abs(Ea_pu) / Xd_new ) * sin(loading_angle) + ( 0.5*abs(Ea_pu)^2 ) * ( (1/Xq_new) - (1/Xd_new) ) * sin(2*loading_angle);
Q_gen_pu = ( abs(Ea_pu) * ( abs(Ei_pu)*cos(loading_angle) - abs(Ea_pu) )/Xd_new ) - ( 0.5*abs(Ea_pu)^2 ) * ( (1/Xq_new) - (1/Xd_new) )*(1-cos(2*loading_angle));
S_gen_pu = P_gen_pu + 1i * Q_gen_pu;

S_transformer_pu = abs(I_pu)^2 * 1i * Xt_pu;
P_transformer_pu = real(S_transformer_pu);
Q_transformer_pu = imag(S_transformer_pu);

S_line_pu = abs(I_pu)^2 * line_impedence_pu;
P_line_pu = real(S_line_pu);
Q_line_pu = imag(S_line_pu);

P_load_pu = Sl_pu * power_factor;
Q_load_pu = sqrt(Sl_pu^2 - P_load_pu^2);
if (lagging == false)
    Q_load_pu = -Q_load_pu;
end

P_consumed_pu = P_transformer_pu + P_line_pu + P_load_pu;
Q_consumed_pu = Q_transformer_pu + Q_line_pu + Q_load_pu;

fprintf("System Information:\n");
% print if load is lagging or leading
if (lagging == true)
    fprintf("Lagging load.\n");
else
    fprintf("Leading load.\n");
end
% print line impedence
fprintf("Line Impedence: %.4f + j %.4f ohms\n", line_resistance_per_km * length, line_reactance_per_km * length);
fprintf("Line Length: %.4f\n", length);

% print voltages (SI) and loading angle
fprintf("\nVOLTAGES AND LOADING ANGLE (System International)\nInternal Voltage: %.4f V\nTerminal Voltage: %.4f V @ %.2f°\nLoading Angle: %.2f°", abs(Ei), abs(Ea), angle(Ea)*180/pi, loading_angle*180/pi);

% print intermediate values (SI)
fprintf("\n\nINTERMEDIATE VALUES (System International)\nLoad Current: %.4f A @ %.2f°\nQuadature Voltage: %.4f V @ %.2f°\nTransformer Voltage: %.4f V @ %.2f°", abs(I2), angle(I2)*180/pi, abs(Eq), angle(Eq)*180/pi, abs(Vt), angle(Vt)*180/pi);

% print voltages (PU) and loading angle
fprintf("\n\nVOLTAGES AND LOADING ANGLE (Per Unit)\nInternal Voltage: %.2f\nTerminal Voltage: %.2f @ %.2f°\nLoading Angle: %.2f°", abs(Ei_pu), abs(Ea_pu), angle(Ea_pu)*180/pi, loading_angle*180/pi);

% print intermediate values (PU)
fprintf("\n\nINTERMEDIATE VALUES (Per Unit)\nLoad Current: %.4f @ %.2f°\nQuadature Voltage: %.4f @ %.2f°\nTransformer Voltage: %.4f @ %.2f°", abs(I_pu), angle(I_pu)*180/pi, abs(Eq_pu), angle(Eq_pu)*180/pi, abs(Vt_pu), angle(Vt_pu)*180/pi);

% print power equilibrium
fprintf("\n\nPOWER EQUILIBIRUM\nGenerated Power = %.4f + j %.4f\nTransformer Power = %.4f + j %.4f\nLine Power = %.4f + j %.4f\nLoad Power = %.4f + j %.4f\nConsumed Active Power = %.4f + %.4f + %.4f = %.4f\nConsumed Reactive Power = %.4f + %.4f + %.4f = %.4f\n", P_gen_pu, Q_gen_pu, P_transformer_pu, Q_transformer_pu, P_line_pu, Q_line_pu, P_load_pu, Q_load_pu, P_transformer_pu, P_line_pu, P_load_pu, P_consumed_pu, Q_transformer_pu, Q_line_pu, Q_load_pu, Q_consumed_pu);