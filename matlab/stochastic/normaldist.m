function n = normaldist(x, mu, sigma)

% the normal distribution pdf n = normaldist(x, mu, sigma)

n = 1./(sigma * sqrt(2 * pi)) .* exp((-(x - mu).^2)./(2.*sigma.^2));