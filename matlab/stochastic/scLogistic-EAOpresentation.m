% scLogistic - simulations with the logistic model


X0 = 0.1;
Xmax = 3;
uMax = 0.5;
alpha = 50;

t = 0:34; % sampled days

tSubSample = 0:0.1:34;

% logistic curve derived from integration
x = X0; 
xSloughing = X0; 
for i = 2:length(tSubSample),
    deltaT = tSubSample(i) - tSubSample(i-1);
    deltaX = deltaT * uMax / Xmax * x(i-1) * ( Xmax - x(i-1) );
    deltaXSloughing = deltaT * uMax / Xmax * xSloughing(i-1) *...
        ( Xmax - xSloughing(i-1) )...
        - rand(1)^alpha*xSloughing(i-1);
    x = [x, x(i-1) + deltaX];
    xSloughing = [xSloughing, xSloughing(i-1) + deltaXSloughing];
end;

figure(2);
set(2, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(tSubSample, x);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);


xObserved = interp1(tSubSample, x, t, 'nearest');

figure(3)
set(3, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(t, xObserved);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);
print('-depsc', 'logistic.eps')

figure(4)
set(4, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(x(1:end-1), x(2:end));
xlabel('x(t)');
ylabel('X(t+{\Delta}t)');
%set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);

figure(5)
set(5, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(xObserved(1:end-1), xObserved(2:end));
xlabel('x(t)');
ylabel('X(t+{\Delta}t)');
print('-depsc', 'logisticPhase.eps')

% with sloughing
figure(6);
set(6, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(tSubSample, xSloughing);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);

xSloughingObserved = interp1(tSubSample, xSloughing, t, 'nearest');

figure(7)
set(7, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(t, xSloughingObserved);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);
print('-depsc', 'sloughingCurve.eps')

figure(8)
set(8, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(xSloughing(1:end-1), xSloughing(2:end));
xlabel('xSloughing(t)');
ylabel('XSloughing(t+{\Delta}t)');
%set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);
print('-depsc', 'sloughingPhaseHighSampling.eps');


figure(9)
set(9, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
%plot(xObserved(1:end-1), xObserved(2:end), 'r--');
% hold on;
plot(xSloughingObserved(1:end-1), xSloughingObserved(2:end));
% hold off;
xlabel('xSloughing(t)');
ylabel('XSloughing(t+{\Delta}t)');
print('-depsc', 'sloughingPhase.eps')


% with sloughing and noise
xSloughingWithNoise = xSloughing + (rand(size(xSloughing))-0.5) * 0.2; 

figure(10)
set(10, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(tSubSample, xSloughingWithNoise);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);

xSloughingWithNoiseObserved =...
    interp1(tSubSample, xSloughingWithNoise, t, 'nearest');

figure(11)
set(11, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(t, xSloughingWithNoiseObserved);
xlabel('t [day]');
ylabel('X');
set(gca, 'XLim', [0 35], 'YLim', [0 3.5]);
print('-depsc', 'sloughingWithNoise.eps')

figure(12)
set(12, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(xSloughingWithNoise(1:end-1), xSloughingWithNoise(2:end));
xlabel('xSloughing(t)');
ylabel('XSloughing(t+{\Delta}t)');
print('-depsc', 'sloughingWithNoisePhaseHighSampling.eps')

figure(13)
set(13, 'Position', [256   477   409   215], 'PaperPositionMode', 'auto');
plot(xSloughingWithNoiseObserved(1:end-1),...
    xSloughingWithNoiseObserved(2:end));
xlabel('xSloughing(t)');
ylabel('XSloughing(t+{\Delta}t)');
print('-depsc', 'sloughingWithNoisePhase.eps')
