L = 0.5; % Length of the enclosure (m)
W = 0.3; % Width of the enclosure (m)
H = 0.4; % Height of the enclosure (m)
Lt = L; Wt = W; % Top surface dimensions
Ti_guess = 90; % Initial guess for internal air temperature (°C)
Tis_guess = 80; % Initial guess for internal wall surface temperature (°C)
Tes_guess = 70; % Initial guess for external surface temperature (°C)
Tamb = 60; % Ambient temperature (°C)
Q = 7; % Heat dissipation (W)
t = 0.01; % Thickness of the enclosure walls (m)
k = 0.13; % Thermal conductivity of the enclosure material (W/m·K)
boltz = 5.67e-8; % Stefan-Boltzmann constant (W/m²·K⁴)
emmis = 0.9; % Emissivity of the enclosure material

% Derived Parameters
A_end = L * W; % Area of the end walls (m²)
L_star = (L * W) / (2 * (W + L)); % Characteristic length for convection (m)
A_vert = 2 * H * (L + W); % Vertical surface area (m²)
A_total = 2 * A_end + A_vert; % Total enclosure surface area (m²)

% Intermediate Calculations for Resistances
H_bi = 0.59 * ((Ti_guess - Tis_guess) / L_star)^0.25;
R_bi = 1 / (H_bi * A_end);

H_be = 0.59 * ((Tes_guess - Tamb) / L_star)^0.25;
R_be = 1 / (H_be * A_end);

H_ti = 1.32 * ((Ti_guess - Tis_guess) / L_star)^0.25;
R_ti = 1 / (H_ti * A_end);

H_te = 1.32 * ((Tes_guess - Tamb) / L_star)^0.25;
R_te = 1 / (H_te * A_end);

H_vi = 1.42 * ((Ti_guess - Tis_guess) / L_star)^0.25;
R_vi = 1 / (H_vi * A_vert);

H_ve = 1.42 * ((Tes_guess - Tamb) / L_star)^0.25;
R_ve = 1 / (H_ve * A_vert);

H_radi = (boltz * (Ti_guess + Tis_guess) * (Ti_guess^2 + Tis_guess^2)) / ((1 / emmis) + ((1 - emmis) / emmis));
R_radi = 1 / (H_radi * A_total);

R_cond = t / (k * A_total);

H_rade = emmis * boltz * (Tes_guess + Tamb) * (Tes_guess^2 + Tamb^2);
R_rade = 1 / (H_rade * A_total);

R_totale = ((1 / R_rade) + (1 / R_be) + (1 / R_te) + (1 / R_ve))^-1;
R_totali = ((1 / R_radi) + (1 / R_bi) + (1 / R_ti) + (1 / R_vi))^-1;

% Additional Thermal Resistance Calculations (from heat sink formulas)
H_hs_conv = 10; % Assumed convection coefficient for heat sink (W/m²·K)
A_hs = 0.04; % Heat sink area (m²)
R_hs_conv = 1 / (H_hs_conv * A_hs);

H_hs_radi = emmis * boltz * ((Ti_guess + Tamb) * (Ti_guess^2 + Tamb^2));
R_hs_radi = 1 / (H_hs_radi * A_hs);

R_hs_total = 1 / (1 / R_hs_conv + 1 / R_hs_radi);

% Define the system of equations
equations = @(T) [
    % Energy balance for Tes (external surface temperature)
    R_totale * Q - (T(3) - Tamb);
    % Energy balance for Tis (internal wall surface temperature)
    R_cond * Q + T(3) - T(2);
    % Energy balance for Ti (internal air temperature)
    R_totali * Q - (T(1) - T(2))];

% Solve the system of equations using fsolve
options = optimset('Display', 'iter');
initial_guesses = [Ti_guess, Tis_guess, Tes_guess];
solution = fsolve(equations, initial_guesses, options);

% Extract the solutions
Ti = solution(1);
Tis = solution(2);
Tes = solution(3);

% Display the results
fprintf('Internal Air Temperature (Ti): %.2f °C\n', Ti);
fprintf('Internal Wall Surface Temperature (Tis): %.2f °C\n', Tis);
fprintf('External Surface Temperature (Tes): %.2f °C\n', Tes);

% Calculate heat source temperature
T_source = Ti + Q * R_hs_total;
fprintf('Heat Source Temperature (T_source): %.2f °C\n', T_source);
