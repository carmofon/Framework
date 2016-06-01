% scTesteOdeNConversion

%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
[t,y] = ode45(@nConversions,[0 3],[0.035 0 0.035]);

%
figure(1);
plot(t,y(:,1),'-',t,y(:,2),'-',t,y(:,3),'-');
set(gca, 'YLim', [0 0.2]);