%{

05.26 - values for s are negative because Alpha < Beta. Either the
parameter values do not numerically make sense or I wrote the equations 
wrong.

Inputs:   fs = vector of f_hat values
          rs = vector of r_hat values

%}

function [s_hats, dfdts, drdts] = HSSrfQuiver(fs, r_hats)

% Model parameter values.
Rs = 0.01;                  % 1/h
rho = 3000;                 % g/L
Ds = 5.76e-6;               % um^2/h
h = 0.05 * 2000;            % um
Q = Ds / h^2;
S0 = 1e-6;                  % g/L
Ys = 0.45;                  % dimensionless
muMax = 0.1;                % 1/h
c = 0.3;                    % dimensionless
Ks = 3.5e-5;                % g/L
k = Ks / S0;
Rdet = 1e-5;                % g/(um^3 * h)
lambda = 0.03;              % 1/h

% Calculate dimensionless lump parameters.
Alpha = (Rs * rho) / (2 * Q * S0);
Beta = (rho * muMax) / (2 * Q * S0 * Ys);
Gamma = (Rdet * h) / (muMax * rho);
Delta = (c * Rs) / muMax;
Epsilon = lambda / muMax;

% Calculate normalized signal concentration for each f hat - r hat pair.
dfdts = zeros(length(r_hats), length(fs));
drdts = dfdts;
s_hats = dfdts;


% no need to use nested loops. Matlab can calculate using vectorized
% notation. This is mush faster (but it requires more memory)

% meshgrid expands the arrays fs and r_hats into matrices
[fs, r_hats] = meshgrid(fs, r_hats);
% then you can calculate using array multplication (.*), 
% division (./) and power (.^2)
s = (Alpha * fs - Beta) .* r_hats + 1;
dfdts = Delta * fs.^2 + (Epsilon - Delta) * fs + Epsilon;
drdts = - Gamma * r_hats.^2 ...
    + 1/2 .* r_hats .* (s ./ (k + s) ...
    - Delta * fs);


% for i = 1:size(dfdts,1)
%     for j = 1:size(dfdts,2)
%         f = fs(j);
%         r = r_hats(i);
%         s = (Alpha * fs(j) - Beta) * r_hats(i) + 1;
%         s_hats(i,j) = s;
%         dfdts(i,j) = Delta * f^2 + (Epsilon - Delta) * f + Epsilon;
%         drdts(i,j) = - Gamma * r^2 ...
%             + 1/2 * r * (s / (k + s) ...
%             - Delta * f);
%     end
% end

quiver(fs, r_hats, dfdts, drdts, 0.25);