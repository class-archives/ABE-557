%% Fermentation Design Algorithm Part I
% Kathryn Atherton

clear;
clc;

%% Starting Variables and Constants
output_rate = 100; % dry solid output rate [lbs/h] (given)
K_s = 0.25; % Michaelis-Menton constant of reaction [g/L] (given)
S_o = 0.2; % initial substrate concentration [g/g-solution] (given)
u_max = 0.5; % maximum growth rate [1/h] (given)
Yxs = 0.5; % cell mass yield based on limiting nutrient [g cells produced / g substrate used] (given)
sol_density = 995.68; % fermentation broth density [g/L] (Geankoplis)
lb_g = 453.592; % conversion from lbs to g
S_final_perc = 0.05; % percent of substrate left at end of fermentation (Assumption 1)
X_o_perc = 0.1; % percent of final concentration of yeast cells at beginning of fermentation (Assumption 2)
fill_ferment_ratio = 0.25; % ratio of fill time to ferment time (Assumption 4)
l_m3 = 1/1000; % conversion from L to m^3
da_dt_ratio = 0.5; % ratio of the agitator to tank diameter (Assumption 7)
working_v_perc = 0.8; % percent of tank that the fermentation broth fills (Assumption 8)
ht_dt_ratio = 3; % ratio of tank diameter to tank height (Assumption 6)

%% Finding Fermentation Time
% Convert the initial substrate concentration from g/g to g/L
S_i = S_o * sol_density; % units = g/g-solution * g-solution/L = g/L
% Find the final substrate concentration
S_f = S_i * S_final_perc; % units = g/L
% Find the final yeast concentration (Equation 1)
X_f = (Yxs * (S_i - S_f)) / (1 - X_o_perc); % units = g/L
% Find the initial yeast concentration 
X_i = X_f * X_o_perc; % units = g/L

% Calculate fermentation time (Equation 2)
ferment = (((K_s * Yxs + S_i * Yxs + X_i) / (Yxs * S_i + X_i)) * log(X_f / X_i) - ((K_s * Yxs) / (Yxs * S_i + X_i)) * log((Yxs * S_i + X_i - X_f) * Yxs * S_i))/u_max; % units = h
fprintf('The fermentation time of the tank is %.2f hours.\n',ferment);
% Calculate the filling time (Assumption 4)
fill = fill_ferment_ratio * ferment; % units = h
% Calculate the emptying time (Assumption 3)
empty = ferment + fill; % units = h

% Calculate volume of one fermentation tank 
volume = (lb_g * output_rate * empty)/X_f; % units = L
fprintf('The volume of one batch is %.2f L.\n',volume);

%% Find Tank Dimensions
% Convert volume from L to m^3 and incoporate working volume (Assumption 8)
t_volume_m = volume * l_m3 / working_v_perc; % units = m^3
t_volume_l = t_volume_m / l_m3; % units = L
% Find diameter of tank (Equation 3, Assumption 5, Assumption 6)
dt = ((t_volume_m * 4) / (pi  * 3))^(1/3); % units = m
% Find agitator diameter (Assumption 7)
da = dt * da_dt_ratio; % units = m
% Find height of tank (Assumption 6)
ht = dt * ht_dt_ratio; % units = m
fprintf('The volume of the tank is %.2f L. The tank has a diameter of %.2f m and a height of %.2f m. The agitator diameter is %.2f m.\n',t_volume_l, dt, ht, da)