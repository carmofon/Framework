% scPlotSolutions

tubeReactor.Y_SX = 0.5; %yield [kgCOD-S/kgCOD-X]
tubeReactor.rho_X = 2e2; %yield [kgCOD-S/m3]
tubeReactor.kDet = 0.5e-15;
tubeReactor.Qv = 14.1145; %volumetric flow [m^3/day]
tubeReactor.eta = 86.4; %dynamic viscosity of water [kg/(m.day)]
tubeReactor.L = 1; %distance from tubeReactor entrance [m]
tubeReactor.Ds = 1.992e-4; %diffusivity of oxygen in water [m^2/day]
tubeReactor.mu_max = 1; %maximum specific growth rate [day^-1]
tubeReactor.R = 0.013; %radius of tubeReactor [m]
tubeReactor.C_S_bulk = 8e-3; %bulk concentration of oxygen [kgO/m^3]


kdet = linspace(1e-17, 1.2e-15, 100); 

for i = 1:length(kdet),
    tubeReactor.kDet = kdet(i);
    [r_f(i), r_delta(i), C_S_surf(i)] =...
        solveSteadyState(tubeReactor);
    tubeReactor.kDet = kdet(i);
    [r_f_d(i), r_delta_d(i), C_S_surf_d(i)] =...
        solveSteadyStateDimensionless3(tubeReactor);
end;

figure(1);
subplot(1,2,1);
plot(kdet, (tubeReactor.R - r_f)*1e6, 'b-', kdet, (tubeReactor.R - r_f_d)*1e6, 'r-');
xlabel('k_{det}');
ylabel('L_f [{\mu}m]');
subplot(1,2,2);
plot(kdet, (r_delta - r_f)*1e6, 'b-', kdet, (r_delta_d - r_f_d)*1e6, 'r-');
xlabel('k_{det}');
ylabel('{\delta} [{\mu}m]');

