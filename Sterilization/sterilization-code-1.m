%% Sterilization Homework
% ABE 557 
% Kathryn Atherton
% September 20, 2018

%% Problem 5.2-3: Unsteady-State Heating of a Stirred Tank
volume = 0.0283; % units = m^3
To_water = 288.8; % units = K
T_steam = 377.6; % units = K
U = 1136; % units = W/m^2.K
area = 0.372; % units = m^2
Tf_water = 338.7; % units = K
cp_water = 4.179; % units = J/g.K
rho_water = 997; % units = kg/m^3
rho_water = rho_water * 1000; % units = kg / m^3 * g / kg = g / m^3
t = (-1 / (U * area)) * cp_water * rho_water * volume * log(abs((T_steam - Tf_water)/(To_water - T_steam))); % units = s
t = t / 3600; % units = s * h / s = h
fprintf('Problem 5.2-3: The time required to heat the water from %.1f K to %.1f K is %.2f hours.\n',To_water,Tf_water,t);

%% Problem 5.4-5: Cooling Beef with Convective Resistance
x = 45.7; % units = mm
x = x / 1000; % units = mm * m / mm
Ti_meat = 37.78 + 273.17; % units = K
T_air = -1.11+ 273.15; % units = K
h = 38.0; % W/m^2.K
k = 0.498; % W/m.K
alpha = 4.464e-4; % m^2/h
m = 4.0; 
slices = 5; 
del_x = x / slices; % units = mm
t = 0.27; % units = h
t_T1 = Ti_meat;
t_T2 = Ti_meat;
t_T3 = Ti_meat;
t_T4 = Ti_meat;
t_T5 = Ti_meat;
ti = 0;
del_t = (del_x ^ 2) / (m * alpha); % units = h
N = h * del_x / k;
while ti < t
    ti = ti + del_t;
    tdelt_T1 = (1/m) * (2 * N * T_air + (m - (2 * N + 2)) * t_T1 + 2 * t_T2); 
    tdelt_T2 = (1/m) * (t_T3 + (m - 2) * t_T2 + t_T1);
    tdelt_T3 = (1/m) * (t_T4 + (m - 2) * t_T3 + t_T2);
    tdelt_T4 = (1/m) * (t_T5 + (m - 2) * t_T4 + t_T3);
    tdelt_T5 = (1/m) * ((m - 2) * t_T5 + 2 * t_T4); 
    t_T1 = tdelt_T1;
    t_T2 = tdelt_T2;
    t_T3 = tdelt_T3;
    t_T4 = tdelt_T4;
    t_T5 = tdelt_T5;
end
x = [0, 1 * del_x * 1000, 2 * del_x * 1000, 3 * del_x * 1000, 4 * del_x * 1000, 5 * del_x * 1000]; % units = mm
T = [T_air - 273.15, t_T1 - 273.15, t_T2 - 273.15, t_T3 - 273.15, t_T4 - 273.15, t_T5 - 273.15]; % units = deg C
plot(x,T);
ylim([T_air - 273.15,Ti_meat - 273.15]);
xlim([0, del_x * 1000]);
xlabel('thickness [mm]')
ylabel('temperature [deg C]')
title('Problem 5.4-5: Temperature vs. Thickness')
fprintf('The temperature at node 1 is %.2f, at node 2 is %.2f, at node 3 is %.2f, at node 4 is %.2f, and at node 5 is %.2 f.\n',T(2), T(3), T(4), T(5), T(6));

%% Problem 9.12-5: Process Time and Numerical Integration
t = [0, 20, 40, 60, 80, 90, 100]; % units = min
T_c = [43.3, 79.3, 96.1, 108.9, 111.1, 107.2, 71.1]; % units = deg C
T_f = [110, 165, 205, 228, 232, 225, 160]; % units = deg F
F0 = 2.60; % units = min
z_c = 10; % units = deg C
z_f = 18; % units = deg F
F0_f = ((10^((T_f(1) - 250)/z_f) + 10^((T_f(2) - 250)/z_f)/2) * (t(2)-t(1))) + ((10^((T_f(2) - 250)/z_f) + 10^((T_f(3) - 250)/z_f)/2) * (t(3)-t(2))) + ((10^((T_f(3) - 250)/z_f) + 10^((T_f(4) - 250)/z_f)/2) * (t(4)-t(3))) + ((10^((T_f(4) - 250)/z_f) + 10^((T_f(5) - 250)/z_f)/2) * (t(5)-t(4))) + ((10^((T_f(5) - 250)/z_f) + 10^((T_f(6) - 250)/z_f)/2) * (t(6)-t(5))) + ((10^((T_f(6) - 250)/z_f) + 10^((T_f(7) - 250)/z_f)/2) * (t(7)-t(6))); % Trapezoidal Method
F0_c = ((10^((T_c(1) - 121.1)/z_c) + 10^((T_c(2) - 121.1)/z_c)/2) * (t(2)-t(1))) + ((10^((T_c(2) - 121.1)/z_c) + 10^((T_c(3) - 121.1)/z_c)/2) * (t(3)-t(2))) + ((10^((T_c(3) - 121.1)/z_c) + 10^((T_c(4) - 121.1)/z_c)/2) * (t(4)-t(3))) + ((10^((T_c(4) - 121.1)/z_c) + 10^((T_c(5) - 121.1)/z_c)/2) * (t(5)-t(4))) + ((10^((T_c(5) - 121.1)/z_c) + 10^((T_c(6) - 121.1)/z_c)/2) * (t(6)-t(5))) + ((10^((T_c(6) - 121.1)/z_c) + 10^((T_c(7) -121.1)/z_c)/2) * (t(7)-t(6))); % Trapezoidal Method
fprintf('Problem 9.12-5: The F0 value in English units is %.2f mins and in SI units is %.2f mins.', F0_f, F0_c);
if ((F0_f + F0_c)/2) >= F0
    fprintf('Therefore, the thermal process is adequate.\n')
else
    fprintf('Therefore the thermal process is inadequate.\n')
end

%% Problem 9.12-6: Sterility Level of Fermentation Medium
t = [0, 10, 20, 25, 30, 35]; % units = min
T_c = [100, 110, 120, 120, 110, 100]; % units = deg C
T_k = T_c + 273.15; % units = K
k = 7.94e38 * exp((-68.7e3) ./ (1.987 .* T)); % units = min^-1
N0 = 1e12; % units = spores
N_final = N0 * exp(-k(6) * t(6)); % units = spores
V = log(N0/N_final);
fprintf("Problem 9.12-6: The final sterility level is %.2d and the V is %.2f.",N_final,V)

%% Problem 9.12-7: Time for Pasteurization of Milk
T = 62.8; % units = deg C
F0 = 9.0; % units = minutes
T1 = 65.6; % units = deg C
z = 5; % units = deg C
t = F0 / (10 ^ ((T-T1)/z)); % units = minutes
fprintf('Problem 9.12-7: The time to pasteurize is %.2f min.\n',t)