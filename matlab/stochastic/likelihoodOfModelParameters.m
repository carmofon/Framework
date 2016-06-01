function [L, params] = likelihoodOfModelParameters(x2, x1, params)

% determine the conditional density value for the array of
% values in x2 given the values in x1


disp('Computing likelihood of parameters in experimental set');
disp(sprintf('mu = %0.2f \t b = %0.2f  \t alpha = %d',...
    params.mu, params.b, params.alpha));

L = 0; % initialize
for i = 1:length(x1),
    l = likelihood(x2(i), x1(i), params);
    if (l == 0)
        L = NaN;
        break;
    else
        L = L + log(l);
    end;
end;

%disp(sprintf('%d of %d (l = %f)', i, length(x1), l));
disp('Finished determining likelihood');
disp(sprintf('L = %f', L));
return;
