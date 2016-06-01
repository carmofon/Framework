% scSolution - analytical solution for the steady state of
% tube reactor

% parameteres
% Y_SX = 1.5; %yield [kgCOD-S/kgCOD-X]
% rho_X = 2.5e5; %yield [kgCOD-S/m3]
% kDet = 1e-15;
% Qv = 3.53; %volumetric flow [m^3/day]
% eta = 86.4; %dynamic viscosity of water [kg/(m.day)]
% L = 1; %distance from tube entrance [m]
% Ds = 1.992e-4; %diffusivity of oxygen in water [m^2/day]
% mu_max = 11.28; %maximum specific growth rate [day^-1]
% R = 0.013; %radius of tube [m]
% C_S_bulk = 8e-3; %bulk concentration of oxygen [kgO/m^3]
Y_SX = 1.5; %yield [kgCOD-S/kgCOD-X]
rho_X = 2.5e2; %yield [kgCOD-S/m3]
kDet = 0.5e-15;
Qv = 3.53; %volumetric flow [m^3/day]
eta = 86.4; %dynamic viscosity of water [kg/(m.day)]
L = 1; %distance from tube entrance [m]
Ds = 1.992e-4; %diffusivity of oxygen in water [m^2/day]
mu_max = 5.28; %maximum specific growth rate [day^-1]
R = 0.013; %radius of tube [m]
C_S_bulk = 4e-3; %bulk concentration of oxygen [kgO/m^3]


N = 10000;
R_f = logspace(log10(R*1e-4),log10(R),N);

R_delta_fun = @(R_f)(sqrt(R_f.^2 + 8*kDet*Qv*eta/mu_max./R_f.^2));
R_delta = R_delta_fun(R_f);


validIndexes = find((R_f < R_delta) & (R_delta < R));
R_f = R_f(validIndexes);
R_delta = R_delta(validIndexes);



%build function of first solution of C at surf
C_S_surf_1 =...
    @(R_f, R_delta)...
    (0.5*Y_SX*mu_max*rho_X/Ds*(R_delta.^2.*(1-log(R_delta./R_f)) - R_f.^2));

C_S_surf_1_vals = C_S_surf_1(R_f, R_delta);
valid_1 = (R_f < R_delta) & (C_S_surf_1_vals > 0) & (C_S_surf_1_vals < 2*C_S_bulk);
C_S_surf_1_vals(~valid_1) = NaN;


%build function of second solution of C at surf
C_S_surf_2 =...
    @(R_f, R_delta)...
    C_S_bulk...
    - sqrt(pi*L/4/Qv/Ds) * 0.5 * Y_SX * mu_max * rho_X *...
    (R_delta.^2 - R_f.^2);

C_S_surf_2_vals = C_S_surf_2(R_f, R_delta);
valid_2 = (R_f < R_delta) & (C_S_surf_2_vals > 0) & (C_S_surf_2_vals < C_S_bulk);
C_S_surf_2_vals(~valid_2) = NaN;


% build solution
f_to_solve = @(R_f, R_delta)...
    (C_S_surf_2(R_f, R_delta) - C_S_surf_1(R_f, R_delta));
f_to_solve_R_f = @(R_f)(f_to_solve(R_f, R_delta_fun(R_f)));

solution = fzero(f_to_solve_R_f, R);


figure(2);
subplot(1,3,1);
plot(R_f, R_delta, 'b-');
xlabel('R_f [m]');
ylabel('R_delta [m]');

subplot(1,3,2);
plot(R_f, C_S_surf_1_vals, 'b-', R_f, C_S_surf_2_vals, 'r-');
xlabel('R_f [m]');
ylabel('C_S^{surf} solution 1 [kgS/m^3]');

subplot(1,3,3);
plot(R_f, f_to_solve_R_f(R_f), 'b-');
xlabel('R_f [m]');
ylabel('f to solve');

