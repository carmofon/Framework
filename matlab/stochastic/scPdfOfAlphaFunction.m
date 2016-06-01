% scPdfOfAlphaFunction

% Determine the probability density of the the alpha function

ALPHA   = 100;
NPOINTS = 10000; 

% generate values using mote-carlo
x = (1 - rand(1, NPOINTS)).^ALPHA;

figure(1);
plot(sort(x));
title('Results for Monte-Carlo sims');

%Estimate PDF
[pdf, xf] = ksdensity(x, 'npoints', 1000, 'kernel', 'normal');

figure(2);
semilogy(xf, pdf);
xlabel('x_f');
ylabel('P.D.F.');
