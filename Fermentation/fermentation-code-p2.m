%% Fermentation Design Algorithm Part II
% Kathryn Atherton
clc;
clear;
%% PART I
%% Finding Fermentation Time
% Starting values, constants, & conversions
S_o = 0.2; % initial substrate concentration [g/g-solution] (given)
sol_density = 995.68; % fermentation broth density [g/L] (Geankoplis)
S_final_perc = 0.05; % percent of substrate left at end of fermentation (Assumption 1)
Yxs = 0.5; % cell mass yield based on limiting nutrient [g cells produced / g substrate used] (given)
X_o_perc = 0.1; % percent of final concentration of yeast cells at beginning of fermentation (Assumption 2)
K_s = 0.25; % Michaelis-Menton constant of reaction [g/L] (given)
u_max = 0.5; % maximum growth rate [1/h] (given)
fill_ferment_ratio = 0.25; % ratio of fill time to ferment time (Assumption 3)
output_rate = 100; % dry solid output rate [lbs/h] (given)
lb_g = 453.592; % conversion from lbs to g [g/lbs] (conversion)

% Convert the initial substrate concentration from g/g to g/L
S_i = S_o * sol_density; % [g/L]
% Find the final substrate concentration
S_f = S_i * S_final_perc; % [g/L]
% Find the final yeast concentration (Equation 1)
X_f = (Yxs * (S_i - S_f)) / (1 - X_o_perc); % [g/L]
% Find the initial yeast concentration 
X_i = X_f * X_o_perc; % [g/L]

% Calculate fermentation time (Equation 2)
ferment = (((K_s * Yxs + S_i * Yxs + X_i) / (Yxs * S_i + X_i)) * log(X_f / X_i) - ((K_s * Yxs) / (Yxs * S_i + X_i)) * log((Yxs * S_i + X_i - X_f) * Yxs * S_i)) / u_max; % [h]
fprintf('The fermentation time of the tank is %.2f hours.\n',ferment);
% Calculate the filling time (Assumption 3)
fill = fill_ferment_ratio * ferment; % [h]
% Calculate the emptying time (Assumption 4)
empty = ferment + fill; % [h]

% Calculate volume of one fermentation tank 
volume = (lb_g * output_rate * empty) / X_f; % [L]
fprintf('The volume of one batch is %.2f L.\n',volume)

%% Find Tank Dimensions
% Starting values, constants, & conversions
l_m3 = 1 / 1000; % conversion from L to m^3 [m^3/L] (conversion)
working_v_perc = 0.8; % percent of tank that the fermentation broth fills (Assumption 5)
da_dt_ratio = 0.33; % ratio of the agitator to tank diameter (Assumption 6)
ht_dt_ratio = 3; % ratio of tank diameter to tank height (Assumption 7)
h_to_s = 3600; % conversion from hours to seconds [s/h] (conversion)

% Convert volume from L to m^3 and incoporate working volume (Assumption 5)
t_volume_m = volume * l_m3 / working_v_perc; % [m^3]
t_volume_l = t_volume_m / l_m3; % [L]
% Find diameter of tank (Equation 3, Assumption 7, Assumption 8)
dt = ((t_volume_m * 4) / (pi  * 3)) ^ (1 / 3); % [m]
% Find agitator diameter (Assumption 6)
da = dt * da_dt_ratio; % [m]
% Find height of tank (Assumption 7)
ht = dt * ht_dt_ratio; % [m]
fprintf('The volume of the tank is %.2f L. The tank has a diameter of %.2f m and a height of %.2f m. The agitator diameter is %.2f m.\n',t_volume_l, dt, ht, da)

% Calculate fill flow rate:
fill_flow = volume / fill; % [L/h]
% Convert to L/s
fill_flow = fill_flow / h_to_s; % [L/s]
fprintf('The fill flow rate is %.2f L/s.\n',fill_flow);

%% Graphing Relationship Between Yeast and Substrate Concentrations Over Time
% Vector length of 100 ensures substrate concentration reaches zero.
% Yeast concentration vectors from 0-200 to find time to reach these concentrations
X_t = zeros(1,100);
for n = 1:length(X_t) 
    X_t(n) = X_i + (n-1);
end
% Fermentation times required to reach yeast concentrations
t_ferm = zeros(1,100);
for n = 1:length(t_ferm)
    t_ferm(n) = (((K_s * Yxs + S_i *Yxs + X_i)/(Yxs * S_i + X_i))*log(X_t(n)/X_i) - ((K_s*Yxs)/(Yxs*S_i + X_i))*log((Yxs*S_i + X_i - X_t(n))/Yxs*S_i))/u_max; % Equation 2)
end
% Substrate concentrations in fermentor for given yeast concentration 
substrate_t = zeros(1,100);
for n = 1:length(substrate_t)
    substrate_t(n) = S_i + ((X_i - X_t(n))/Yxs); % (Equation 1)
end
% Plots yeast and substrate concentration v time
figure
hold on
plot(t_ferm, X_t)
plot(t_ferm, substrate_t)
plot(zeros(1,100)+t_ferm(95),substrate_t, 'r')
xlabel('time (hours)')
ylabel('Concentration (g/L)')
legend('Yeast Concentration', 'Substrate Concentration', 'Location', 'West') 
%% PART II
%% Solving for kLa
% Starting values, constants, & conversions
p_out = 0.1; % partial pressure of oxygen leaving tank [atm O2/atm total] (Assumption 9)
h = 4.75e-4; % Henry's constant for oxygen in water at 30C [atm/mol fraction] (Geankoplis)
m_water = 18.015; % molar mass of water [g/mol] (constant)
qo2 = 8; % respiration rate [mmol O2/g-dw-hr] (given)
c_crit = 1.5 / 32; % lowest oxygen concentration needed by yeast [mmol O2/L] (given)

% calculate c*: soluble oxygen concentration in fermentation broth (from p_out)
c_star = (p_out / h) * (1 / m_water) * sol_density * 1000; % [mmol O2/L]
% solve for kLa: mass transfer coefficient of O2 from air bubbles to fermentation broth (Equation 4)
kLa = (X_f * qo2) / (c_star - c_crit); % [1/h]
fprintf('The mass transfer coefficient (kLa) is %.6f h^-1.\n', kLa)

%% Finding Power Requirement for Agitator
% Starting values, constants, & conversions
n_re = 100000; % Reynolds number (Assumption 10)
mu_water = 7.98e-4; % viscosity of fermentation broth [Pa.s] (Assumption 11)
s_to_min = 60; % conversion of seconds to minutes [s/min] (conversion)
p_num = 5.9; % Power number (Assumption 10)
r = 0.08206; % universal gas constant [L-atm/K-mol] (constant)
t = 30 + 273.15; % temperature of fermentation tank [K] (given)
p_in = 0.21; % partial pressure of oxygen entering tank [atm O2/atm total] (constant in air)

% find rotational speed of impeller in rps (Equation 5)
n_s = (n_re * mu_water) / (sol_density * (da ^ 2)); % [rps]
n_m = n_s * s_to_min; % [rpm]
fprintf('The impeller rotational speed is %0.2f rpm.\n',n_m);
% solve for volumetric flow rate of oxygen
% power of ungassed system (Equation 6)
pu = p_num * sol_density * da^5 * n_s^3; 
% cross-sectional area of fermentation tank (Assumption 8)
a = pi * (dt/2)^2; % [m^2]
% volumetric flow rate of oxygen from ideal gas law and mass balance
q = (((X_f * qo2) * r * t) / (p_in - p_out)) * (volume * l_m3) * (1/1000) * (1/h_to_s); % [m^3/s]
fprintf('The aeration rate is %0.2f m^3/s.\n', q);
% solve for power of gassed system (Equation 7)
pg = 0.354 * (((pu ^ 2) * n_m * (da ^ 3)) / (q / ((volume * l_m3) ^ 0.56))) ^ 0.45;
fprintf('The power requirement for the agitator is %.2f W.\n',pg);

%% Heat Exchanger
% Starting values, constants, & conversions
cp = 4.184; % specific heat of water [J/K] (constant T = 30C)
t_out = 25 + 273.15; % water temperature exiting heat exchanger [K] (Assumption 12)
t_in = 15 + 273.15; % water temperature entering heat exchanger [K] (Assumption 12
velocity = 1.7; % velocity of water flowing through heat exchanger pipe [m/s] (Assumption 13)
wall_thickness = 8.18/1000; % pipe wall thickness [m] (Geankoplis)
k_water = 0.6; % thermal conductivity of water [W/m-K] (Geankoplis)
k_steel = 16; % thermal condutivity of stainless steel at 20C [W/m-K] (Geankoplis, Assumption 14)

% calculate heat generated by cells
q_heat = 0.12 * (qo2 * X_f);  % [kcal/L-h] (Equation 8) 
% convert to kJ/s
q_kj = q_heat * volume * cp / h_to_s; % [kJ/s]
fprintf('The heat generated from oxygen consumption is %.2f kJ/s.\n',q_kj);
% calculate mass flow rate of water within pipes of heat exchanger
m_w = (q_kj)/((t_out - t_in) * cp); % [kg/s]
% convert mass flow rate to cross sectional area of pipe (Assumption 15)
pipe_area = m_w * (1 / sol_density) * (1 / velocity); % [m^2]
% find inner pipe diameter from area
dp_in = 2 * (pipe_area / pi) ^ 0.5; % [m]
% find outer diameter
dp_out = dp_in + 2 * wall_thickness; % [m]
% heat transfer coefficient between pipe and broth 
npr = cp * mu_water * 1000 / k_water; % (Equation 9)
ho = (0.027 * k_water * (da * n_s * sol_density / mu_water) ^ 0.8 * npr ^ (1/3)) / dp_out; % [W/m^2-K] (Equation 10)
hi = (0.027 * k_water * n_re ^ 0.8 * npr ^ (1/3)) / dp_in; % [W/m^2-K] (Equation 10)
% calculate heat transfer for heat exchanger (Equation 11)
u = 1/((1/ho) + (wall_thickness/k_steel) + (1/hi)); % []
% calculate temperature drop of heat exchanger
delt = ((t - t_in) + (t - t_out)) / 2; % [K]
% calculate surface area of heat exchanger to control temperature of fermentation broth in tank (Equation 12)
a_heat_exchanger = (q_kj * 1000) / (u * delt); % [m^2]
fprintf('Surface area of heat exchanger is %.2f m^2.\n',a_heat_exchanger);
l_heat_exchanger = a_heat_exchanger / (pi * dp_out); % [m]

%% Determine Air Flow Rate
% Starting values, constants, & conversions
v_air = 16; % velocity of air through the pipe [m/s] (Assumption 16)
air_pipe_length = pi * 1.5 * dt; % length of the air pipe [m] (Assumption 17)
globe_valve_ld = 300; % (Geankoplis, Table 2.10-1)
roughness = 4.6e-5; % (Geankoplis, Figure 2.10-3)
air_density = 1.177; % air density [g/L] (constant at T = 30)
air_viscosity = 1.846e-5; % viscosity of air [Pa.s] (constant at T = 30)
g = 9.81; % gravitational constant [m/s^2] (constant)
k_sparger = 0.5; % constant for sparger equation (Green & Perry)
j_l_atm = 101.33; % conversion from J to L-atm [J/L-atm] (conversion)

% calculate pipe area
air_pipe_area = (q * volume * l_m3) / v_air; % [m^2] (Assumption 15)
% calculate diameter of air pipe
air_pipe_diameter = ((4 / pi) * air_pipe_area) ^ 0.5; % [m] (Assumption 15)
fprintf('The diameter of the air pipe is %.2f m.\n',air_pipe_diameter);
% calculate pipe equivalent length
pipe_ld = air_pipe_length/air_pipe_diameter;
% find total pipe equivalent length
total_pipe_ld = globe_valve_ld + pipe_ld;
% find equivalent roughness
roughness_d = roughness/air_pipe_diameter;
% calculate air Reynolds number
n_re_air = air_pipe_diameter * air_density * v_air / air_viscosity; % (Equation 5)
fanning = 16 / n_re_air; % (Equation 13)
% find friction from entering pipe (Equation 14)
f_enter = 0.55 * ((v_air ^ 2) / 2); % [J/kg]
fprintf('The friction from entering the pipe is %.2f J/kg.\n', f_enter);
% find friction from pipe
f_pipe_pump = 4 * fanning * air_density * total_pipe_ld * ((v_air ^ 2) / 2); % [J/kg] (Equation 15)
fprintf('The friction from the pipe is %.2f J/kg.\n',f_pipe_pump);
% find friction from sparger (Equation 16)
f_sparger = ((4 * fanning * dt / (3 * air_pipe_diameter)) * 2 * k_sparger) * ((air_density * v_air ^ 2) / 2); % [J/kg] (Equation 21)
fprintf('The friction from the sparger is %.2f J/kg.\n', f_sparger);
% find friction from exit (Equation 17)
f_exit = ((v_air ^ 2) / 2); % [J/kg]
fprintf('The friction from exiting the pipe is %.2f J/kg.\n', f_exit);
% find total friction
f_total = f_enter + f_pipe_pump + f_sparger + f_exit; 
fprintf('The total friction is %.2f J/kg.\n',f_total);
% convert friction to pressure drop
p_drop_total = f_total / j_l_atm * air_density /1000; % [atm]

%% Pump Power Requirements
% Starting values, constants, & conversions
cp_air = 1.005; % heat capacity of air with constant pressure [kJ/kg-K] (Raggett)
cv_air = 0.718; % heat capacity of air with constant volume [kJ/kg-K] (Raggett)
mw_air = 29/1000; % molecular weight of air [g/mol] (constant)
r = 8.314; % universal gas constant [J/mol-K] (constant)
p_in_pump = 1; % suction pressure of pump [atm] (Assumption 18)
pa_atm = 9.86923e-6; % conversion from Pa to atm [atm/Pa] (conversion)
atm_kpa = 101.325; % conversion from atm to kPa [kPa / atm] (conversion)

% Calculate discharge pressure of pump
p_out_pump = p_in_pump + p_drop_total; % [atm]
% convert to kPa
p_kpa = p_out_pump * atm_kpa; % [kPa]
fprintf('The pump must output %.2f kPa of pressure.\n',p_kpa);
% calculate power requirement (Equation 18)
power_pump = (((p_out_pump - p_in_pump) / air_density) / (ferment * h_to_s)) + (f_total * air_density / l_m3 * air_pipe_area * v_air / 1000); % [J/s]
fprintf('The power requirement for the pump is %.2f W.\n', power_pump);

%% Effect of Agitation Speed on System
n_changing = [0:0.01:2]; % range of agitation speeds [rps] (assuming Re is constant for calculation purposes)
pu_changing = p_num .* sol_density .* da^5 .* (n_changing .^3); % power ungassed system [W]
pg = 0.354 * (((pu_changing .^ 2) .* (n_changing * s_to_min) .* (da ^ 3)) ./ (q / ((volume * l_m3) ^ 0.56))) .^ 0.45; % power gassed system [W]
% Plot power vs n
figure
hold on
plot(n_changing, pg)
xlabel('Agitation Speed [rps]')
ylabel('Power of Gassed System [W]')

%% Effect of Air Flow Rate on System
v_air_changing = [9:0.5:36]; % range of air velocity from Geankoplis Table 2.10-3 [m/s]
air_pipe_area_changing = (q * volume * l_m3) ./ v_air_changing; % [m^2] (Assumption 15)
air_pipe_diameter_changing = ((4 / pi) .* air_pipe_area) .^ 0.5; % [m] (Assumption 15)
pipe_ld_changing = air_pipe_length ./ air_pipe_diameter_changing;
total_pipe_ld_changing = globe_valve_ld + pipe_ld_changing;
roughness_d_changing = roughness ./ air_pipe_diameter_changing;
n_re_air_changing = air_pipe_diameter_changing .* air_density .* v_air_changing ./ air_viscosity; % (Equation 5)
fanning_changing = 16 ./ n_re_air_changing; % (Equation 13)
f_enter_changing = 0.55 .* ((v_air_changing .^ 2) ./ 2); % [J/kg]
f_pipe_pump_changing = 4 .* fanning_changing .* air_density .* total_pipe_ld_changing .* ((v_air_changing .^ 2) ./ 2); % [J/kg] (Equation 15)
f_sparger_changing = ((4 .* fanning_changing .* dt ./ (3 .* air_pipe_diameter_changing)) .* 2 .* k_sparger) .* ((air_density .* v_air_changing .^ 2) ./ 2); % [J/kg] (Equation 21)
f_exit_changing = ((v_air_changing .^ 2) ./ 2); % [J/kg]
f_total_changing = f_enter_changing + f_pipe_pump_changing + f_sparger_changing + f_exit_changing; 
p_drop_total_changing = f_total_changing ./ j_l_atm .* air_density ./1000; % [atm]
p_out_pump_changing = p_in_pump + p_drop_total_changing; % [atm]
p_kpa_changing = p_out_pump_changing .* atm_kpa; % [kPa]
power_pump_changing = (((p_out_pump_changing - p_in_pump) ./ air_density) ./ (ferment * h_to_s)) + (f_total_changing .* air_density ./ l_m3 .* air_pipe_area_changing .* v_air_changing ./ 1000); % [J/s]
% Plots power vs air velocity
figure
hold on
plot(v_air_changing, power_pump_changing)
xlabel('Air Velocity [m/s]')
ylabel('Pump Power Requirement [W]')