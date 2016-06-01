%{
Numerical simulation of dr/dt and df/dt [DIMENSIONED] for the HSS
approximation.
%}

function [t, y] = numericalSim(tspan, lambdaX, cX)

% Set initial cooperator fraction and tumor radius
f_ini = 0;
r_ini = 1e-3;

[t, y] = ode45(@derivatives, tspan, [f_ini, r_ini], [], lambdaX, cX);
% [t, y] = ode45(@derivatives, tspan, f_ini, [], lambdaX, cX);

