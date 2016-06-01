% scPlotSolutions

Y_SX = 0.5; %yield [kgCOD-S/kgCOD-X]
rho_X = 2e2; %yield [kgCOD-S/m3]
kDet = 0.5e-15;
Qv = 14.1145; %volumetric flow [m^3/day]
eta = 86.4; %dynamic viscosity of water [kg/(m.day)]
L = 1; %distance from tubeReactor entrance [m]
Ds = 1.992e-4; %diffusivity of oxygen in water [m^2/day]
mu_max = 1; %maximum specific growth rate [day^-1]
R = 0.013; %radius of tubeReactor [m]
C_S_bulk = 8e-3; %bulk concentration of oxygen [kgO/m^3]


%kDet = linspace(1e-17, 1.2e-15, 50);
C_S_bulk = linspace(0, 8e-3, 50);
L = linspace(0, 2, 50);

[L, C_S_bulk] = meshgrid(L, C_S_bulk);

% create the dimensionles groups
c1 = zeros(size(L));
c2 = zeros(size(L));
c3 = zeros(size(L));


c1 = 8.*kDet.*Qv.*eta./mu_max./R.^4;
c2 = 2.*Ds.*C_S_bulk./(Y_SX.*mu_max.*rho_X.*R.^2);
c3 = sqrt(Ds.*pi.*L./Qv);

r_f = zeros(size(L));
r_delta = zeros(size(L));
C_S_surf = zeros(size(L));
solutionType = zeros(size(L));

for i = 1:length(L(:)),
    % solve
    [r_f(i), r_delta(i), C_S_surf(i), solutionType(i)] =...
        solveSteadyStateDimensionless(c1, c2(i), c3(i));
end;

figure(1);
set(1, 'Position', [22   178   995   443]);
subplot(1, 3, 1);
pcolor(c3, c2, r_f)
xlabel('dimensionless group 3');
ylabel('dimensionless group 2');
colorbar;
%axis equal;
subplot(1, 3, 2);
pcolor(c3, c2, solutionType)
xlabel('dimensionless group 3');
ylabel('dimensionless group 2');
colorbar;
%axis equal;
subplot(1, 3, 3);
pcolor(c3, c2, C_S_surf)
xlabel('dimensionless group 3');
ylabel('dimensionless group 2');
colorbar;

nanr_f = zeros(size(r_f));
nanr_f(isnan(r_f)) = 1;

figure(2);
pcolor(c3, c2, nanr_f);
xlabel('dimensionless group 3');
ylabel('dimensionless group 2');
%axis equal;

figure(3);
contour(c3, c2, r_f)
xlabel('dimensionless group 3');
ylabel('dimensionless group 2');
zlabel('R_f');



