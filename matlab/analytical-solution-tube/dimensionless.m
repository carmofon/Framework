%dimensionless

Y_SX = 1.5; %yield [kgCOD-S/kgCOD-X]
rho_X = 2.5e2; %yield [kgCOD-S/m3]
kDet = 0.5e-15;
Qv = 3.53; %volumetric flow [m^3/day]
eta = 86.4; %dynamic viscosity of water [kg/(m.day)]
L = 1; %distance from tubeReactor entrance [m]
Ds = 1.992e-4; %diffusivity of oxygen in water [m^2/day]
mu_max = 5.28; %maximum specific growth rate [day^-1]
R = 0.013; %radius of tubeReactor [m]
C_S_bulk = 4e-3; %bulk concentration of oxygen [kgO/m^3]

c1 = 8*kDet*Qv*eta/mu_max/R^4;
c2 = 2*Ds*C_S_bulk/(Y_SX*mu_max*rho_X*R^2);
c3 = sqrt(Ds*pi*L/Qv);


syms rf rd;
eq1 = rd - sqrt(rf-c1/rf);
eq2 = (rf^2-rf-c1/rf)*(1-c3)+(rf+8*c1/rf)*log(1+8*c1/rf^2)-c2;


f = @(x)((x.^2-x-c1./x).*(1-c3)+(x+8*c1./x).*log(1+8*c1./x.^2)-c2);

%solution = solve(eq2, rf);
solnumeric = fzero(f, 1);

rf_array = linspace(-1,2,10000);
figure(3)
ezplot(eq2);
hold on;
plot(rf_array, f(rf_array), 'r-', solnumeric, f(solnumeric), 'k+');
hold off;
