function dydt = derivatives(t, y, lambdaX, cX)

f = y(1);
r = y(2);

% Model parameter values.
Rs = 0.01;                  % 1/h
rho = 3000;                 % g/L
Ds = 5.76e-6;               % um^2/h
h = 0.05 * 2000;            % um
Q = Ds / h^2;
S0 = 1e-6;                  % g/L
Ys = 0.45;                  % dimensionless
muMax = 0.1;                % 1/h
c = 20;                     % dimensionless
Ks = 3.5e-5;                % g/L
k = Ks / S0;
% Rdet = 1e-5;                % g/(um^3 * h)
Rdet = 5;
% lambda = 0.1;               % 1/h
lambda = 0;

lambda = lambda * lambdaX;
c = c * cX;

% Calculate dimensionless lump parameters.
Alpha = (Rs * rho) / (2 * Q * S0);
Beta = (rho * muMax) / (2 * Q * S0 * Ys);

s = (Alpha * f - Beta) * r/h + 1;

dfdt = (lambda - c*Rs*f) * (1 - f);

% drdt = -1 * Rdet / rho * r^2 ...
%     + muMax * (1/2) * r * s / (k + s) ...
%     - c * Rs * f * r / 2;

% Approximate S0 >> Ks
drdt = -1 * Rdet / rho * r^2 ...
        + (1/2) * muMax * r ...
        - c * Rs * f * r / 2;

dydt = [dfdt; drdt];

% dydt = dfdt;