function x2 = sdeNextStep(x1, params)

% implements the SDE for biofilm development based on the logistic
% equation with stochastic sloughing (power law)

% extract parameters from the params data structure
mu = params.mu;
b = params.b;
a = params.alpha;
n = params.integrationStepsNumber;
dt = params.integrationStepSize;

% generate the random numbers
randomNumber = rand(n, length(x1(:)));
f = (1 - randomNumber).^a;

x2sim = zeros(n, length(x1(:)));

x2sim(1, :) = x1(:);

for i = 2:n,
    x2sim(i, :) = x2sim(i-1, :) + dt * ( mu*x2sim(i-1, :)...
         - b * x2sim(i-1, :).^2) - f(i-1, :) .* x2sim(i-1, :);
end;

x2 = x2sim(end, :);