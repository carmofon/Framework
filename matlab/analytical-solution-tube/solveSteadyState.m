function [r_f, r_delta, C_S_surf] = solveSteadyState(tubeReactor)

Y_SX = tubeReactor.Y_SX; %yield [kgCOD-S/kgCOD-X]
rho_X = tubeReactor.rho_X; %yield [kgCOD-S/m3]
kDet = tubeReactor.kDet;
Qv = tubeReactor.Qv; %volumetric flow [m^3/day]
eta = tubeReactor.eta; %dynamic viscosity of water [kg/(m.day)]
L = tubeReactor.L; %distance from tubeReactor entrance [m]
Ds = tubeReactor.Ds; %diffusivity of oxygen in water [m^2/day]
mu_max = tubeReactor.mu_max; %maximum specific growth rate [day^-1]
R = tubeReactor.R; %radius of tubeReactor [m]
C_S_bulk = tubeReactor.C_S_bulk; %bulk concentration of oxygen [kgO/m^3]

%build the r_delta
R_delta_fun = @(R_f)(sqrt(R_f.^2 + 8*kDet*Qv*eta/mu_max./R_f.^2));


%build function of first solution of C at surf
C_S_surf_1 =...
    @(R_f, R_delta)...
    (0.5*Y_SX*mu_max*rho_X/Ds*(R_delta.^2.*(1-log(R_delta./R_f)) - R_f.^2));

%build function of second solution of C at surf
C_S_surf_2 =...
    @(R_f, R_delta)...
    C_S_bulk...
    - sqrt(pi*L/4/Qv/Ds) * 0.5 * Y_SX * mu_max * rho_X *...
    (R_delta.^2 - R_f.^2);


% build equation to solve
f_to_solve = @(R_f, R_delta)...
    (C_S_surf_2(R_f, R_delta) - C_S_surf_1(R_f, R_delta));
f_to_solve_R_f = @(R_f)(f_to_solve(R_f, R_delta_fun(R_f)));

% compute the solutions
r_f = fzero(f_to_solve_R_f, R);
r_delta = R_delta_fun(r_f);
C_S_surf = C_S_surf_2(r_f, r_delta);

% determine if this is valid solution
if isnan(r_f)
    return
elseif (r_f > R)
    warning ('invalid solution r_f > R');
    r_f = NaN;
    r_delta = NaN;
    C_S_surf = NaN;
elseif (r_f > r_delta)
    warning ('invalid solution r_f > r_delta');
    r_f = NaN;
    r_delta = NaN;
    C_S_surf = NaN;
elseif (r_delta > R)
    warning ('fully penetrated biofilm');
end

