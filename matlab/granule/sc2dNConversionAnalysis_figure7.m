% sc2dNConversionAnalysis_figure7

addpath('../ModelResultBrowser');
basedir = 'e:\jxavier';
resultPath = [basedir '/results'];

dirName = [resultPath '/AnaerobicFeedDO40Precision0_01MaxF0_95_CN50'];

results = getResultsFromDirectory(dirName);
% the size of the system side
results.sizeX = 1600;



%% -------read data from current variable
%iterationsToShow = (64547:64570) + 3;
%iterationsToShow = (75354:75385) + 3;
iterationsToShow = (137671:137693) + 3;

i = iterationsToShow(end - 7);
%for i = iterationsToShow(end-3),
%for i = iterationsToShow(2),
results.iteration.current = i;

% read relevant solutes
results.solutes.show = 1;
results.solids.show = 0;

results.solutes.current = 1;
ammonium = readVariable(results);

results.solutes.current = 3;
nitrite = readVariable(results);

results.solutes.current = 2;
nitrate = readVariable(results);

results.solutes.current = 4;
oxygen = readVariable(results);

results.solutes.current = 6;
acetate = readVariable(results);

results.solutes.current = 5;
phosphate = readVariable(results);

% read relevant solids
results.solids.show = 1;
results.solutes.show = 0;

results.solids.current = 2;
activeNH = readVariable(results);

results.solids.current = 3;
activeNO = readVariable(results);

results.solids.current = 1;
activeH = readVariable(results);

results.solids.current = 4;
activePAO = readVariable(results);

results.solids.current = 10;
phb = readVariable(results);

results.solids.current = 11;
polyp  = readVariable(results);

results.solids.current = 5;
glycogen  = readVariable(results);

% plot the concentrations
%draw2DSubplot(ammonium, 1, i - iterationsToShow(1) + 1);
%draw2DSubplot(nitrite, 2, i - iterationsToShow(1) + 1);
%draw2DSubplot(nitrate, 3, i - iterationsToShow(1) + 1);

% compute the rate of nitrification rNH4
uMaxXNH = 0.0167; % [gCOD-XNH/gCOD-XNH/h] Merle
Y_XNH_NH4 = 0.150; % [gCOD-XNH/gN]
K_XNH_O2 = 0.3e-3; % [gO2/L] Merle
K_XNH_NH4 = 2.4e-3; % [gN/L]

rNH4 = 1/Y_XNH_NH4 * uMaxXNH .* saturation(ammonium, K_XNH_NH4)...
    .* saturation(oxygen, K_XNH_O2) .* activeNH;
%draw2DSubplot(rNH4, 4, i - iterationsToShow(1) + 1);

% compute the rate of convertion of NO2 to NO3
uMaxNO = 0.0458; % [gCOD-XNO/gCOD-XNO/h] Merle
Y_XNO_NO2 = 0.041; % [gCOD-XNO/gN];
K_XNO_O2 = 0.1e-3; % [gO2/L] Merle
K_XNO_NO2 = 0.238e-3; % [gN/L] Merle

r2 = 1/Y_XNO_NO2 * uMaxNO .* saturation(nitrite, K_XNO_NO2)...
    .* saturation(oxygen, K_XNO_O2) .* activeNO;


% compute the rate of denitrification by H rNO3_H
uMaxH = 0.47; % [gCOD-XH/gCOD-XH/h] xavier
K_XH_S = 4e-3; % [gCOD-S/L] xavier (B&B)
Y_XH_S = 0.6645; % [gCOD-XH/gS] xavier (j. biofilms)
K_XH_O2 = 3.5e-4; % [gO2/L] xavier (B&B)
K_XH_NO3 = 1e-3; % [gN/L] same as for PAO
etaNO3 = 0.2; % Merle

rNO3_H = 1/3.43*(1-1/Y_XH_S) * etaNO3 * uMaxH .* saturation(nitrate, K_XH_NO3)...
    .* inhibition(oxygen, K_XH_O2) .* saturation(acetate, K_XH_S) .* activeH;
%draw2DSubplot(rNO3_H, 5, i - iterationsToShow(1) + 1);

%----------------------------------------------------
% compute the rate of denitrification by PAO rNO3_PAO_PHB
YPhbAnoxicNO3PAO = 1.7; % [gCOD-PHB/gCOD-XPAO];
kPhbPAO = 0.23; % [gCOH-PHB/gCOD-PHB/h];
KPAONO3 = 1e-3; % [gN/L]
KPAOO2 = 2e-4; % [gO2/L]
KfPHA_P = 0.33; % [gCOD-PHA/gCOD-XPAO]

rNO3_PAO_PHB = 1/2.86 * (1/YPhbAnoxicNO3PAO-1) * kPhbPAO * etaNO3 * ...
    saturation(nitrate, KPAONO3) .* inhibition(oxygen, KPAOO2) .*...
    saturationFromFraction(phb, activePAO, KfPHA_P) .* activePAO;

rNO3_PAO_PHB(isnan(rNO3_PAO_PHB)) = 0;
%draw2DSubplot(rNO3_PAO_PHB, 6, i - iterationsToShow(1) + 1);

% compute the rate of denitrification by PAO rNO3_PAO_PP
YPpAnoxicNO3PAO = 3.02; % [gP/gCOD] Merle
kPPPAO = 0.0167; % [gP/gCOD-PAO/h];
KPP_P = 0.01; % [gP/gCOD-XPAO];
fPPmaxPAO = 0.65; % [gP/gCOD-XPAO] Merle
KPAOPO4 = 1e-6; % [gP/L];

rNO3_PAO_PP = -1/2.86/YPpAnoxicNO3PAO * kPPPAO * etaNO3 *...
    saturation(nitrate, KPAONO3) .* inhibition(oxygen, KPAOO2) .*...
    inhibitionFromFractionCapacity(polyp, activePAO, KPP_P, fPPmaxPAO) .*...
    saturation(phosphate, KPAOPO4) .* activePAO;

rNO3_PAO_PP(isnan(rNO3_PAO_PP)) = 0;

%draw2DSubplot(rNO3_PAO_PP, 7, i - iterationsToShow(1) + 1);

% compute the rate of denitrification by PAO rNO3_PAO_GLY
YGlyAnoxicNO3PAO = 1.18; % [gCOD-Gly/gCOD] Merle
kGlyPAO = 0.93; % [gCOD-Gly/gCOD-PHB/h] Merle
fGLYmaxPAO = 0.5; % [gCOD-Gly/gCOD-PAO] Merle

rNO3_PAO_GLY = 1/2.86 * (1 - 1/YGlyAnoxicNO3PAO) * kGlyPAO * etaNO3 .*...
    saturation(nitrate, KPAONO3) .* inhibition(oxygen, KPAOO2) .*...
    saturationFromFraction(phb, activePAO, KfPHA_P) .*...
    inhibitionFromFractionCapacity(glycogen, activePAO, 0.01, fGLYmaxPAO).*...
    activePAO;

rNO3_PAO_GLY(isnan(rNO3_PAO_GLY)) = 0;

%draw2DSubplot(rNO3_PAO_GLY, 8, i - iterationsToShow(1) + 1);

%----------------------------------------------------
% compute the rate of denitrification by PAO rNO3_PAO_PHB
YPhbAnoxicNO2PAO = 1.7; % [gCOD-PHB/gCOD-XPAO];
etaNO2 = 0.2; % Merle
KPAONO2 = 1e-3; % [gN/L]

rNO2_PAO_PHB = 1/1.71 * (1/YPhbAnoxicNO3PAO-1) * kPhbPAO * etaNO2 * ...
    saturation(nitrite, KPAONO2) .* inhibition(oxygen, KPAOO2) .*...
    saturationFromFraction(phb, activePAO, KfPHA_P) .* activePAO;

rNO2_PAO_PHB(isnan(rNO2_PAO_PHB)) = 0;
%draw2DSubplot(rNO2_PAO_PHB, 9, i - iterationsToShow(1) + 1);

% compute the rate of denitrification by PAO rNO2_PAO_PP
YPpAnoxicNO2PAO = 3.02; % [gP/gCOD] Merle

rNO2_PAO_PP = -1/1.71/YPpAnoxicNO2PAO * kPPPAO * etaNO2 *...
    saturation(nitrite, KPAONO2) .* inhibition(oxygen, KPAOO2) .*...
    inhibitionFromFractionCapacity(polyp, activePAO, KPP_P, fPPmaxPAO) .*...
    saturation(phosphate, KPAOPO4) .* activePAO;

rNO2_PAO_PP(isnan(rNO2_PAO_PP)) = 0;

%draw2DSubplot(rNO2_PAO_PP, 10, i - iterationsToShow(1) + 1);

% compute the rate of denitrification by PAO rNO2_PAO_GLY
YGlyAnoxicNO2PAO = 1.18; % [gCOD-Gly/gCOD] Merle

rNO2_PAO_GLY = 1/1.71 * (1 - 1/YGlyAnoxicNO2PAO) * kGlyPAO * etaNO2 .*...
    saturation(nitrite, KPAONO2) .* inhibition(oxygen, KPAOO2) .*...
    saturationFromFraction(phb, activePAO, KfPHA_P) .*...
    inhibitionFromFractionCapacity(glycogen, activePAO, 0.01, fGLYmaxPAO).*...
    activePAO;

rNO2_PAO_GLY(isnan(rNO2_PAO_GLY)) = 0;

%draw2DSubplot(rNO2_PAO_GLY, 11, i - iterationsToShow(1) + 1);

%%
imrgb = rNH4;
imrgb(:,:,2) = (-rNO3_PAO_PHB-rNO3_PAO_GLY-rNO2_PAO_PHB-rNO2_PAO_GLY) * 2;
imrgb(:,:,3) = zeros(size(rNH4));

figure(12)
imrgb = uint8((imrgb - min(imrgb(:)))./(max(imrgb(:)) - min(imrgb(:))) * 255);

imshow(imrgb);


%% N-removal

figure(100);
set(100, 'Position', [3         262        1276         662]);

subplot(2,4,1);
contourData2D(ammonium, 'ammonium [gN.L^{-1}]');

subplot(2,4,2);
contourData2D(nitrite, 'nitrite [gN.L^{-1}]');

subplot(2,4,3);
contourData2D(nitrate, 'nitrate [gN.L^{-1}]');

subplot(2,4,4);
contourData2D(oxygen, 'oxygen [gO.L^{-1}]');

r1 = rNH4;
r2 = r2;
r3 = (rNO2_PAO_PHB + rNO2_PAO_GLY);
r4 = (rNO3_PAO_PHB + rNO3_PAO_GLY);

subplot(2,4,5);
rasterData2D(r1, 'rate R1 [gN.L^{-1}h^{-1}]');

subplot(2,4,6);
rasterData2D(r2, 'rate R2 [gN.L^{-1}h^{-1}]');

subplot(2,4,7);
rasterData2D(r3, 'rate R3 [gN.L^{-1}h^{-1}]');

subplot(2,4,8);
rasterData2D(r4, 'rate R4 [gN.L^{-1}h^{-1}]');

%% Gradients of concentration
%- Ammonium
figure(13);
contour(ammonium);
hold on;
[dAmmoniumX, dAmmoniumY] = gradient(ammonium);
quiver(dAmmoniumX, dAmmoniumY, 2);
hold off;
axis equal tight;
%- Nitrite
figure(14);
contour(nitrite);
hold on;
[dNitriteX, dNitriteY] = gradient(nitrite);
quiver(dNitriteX, dNitriteY, 2);
hold off;
axis equal tight;
%- Nitrate
figure(15);
contour(nitrate);
hold on;
[dNitrateX, dNitrateY] = gradient(nitrate);
quiver(dNitrateX, dNitrateY, 2);
hold off;
axis equal tight;

%% draw the particles in grayscale
%figure(16);
%drawPartilcesInGrayscale (results);


%end;
tdata = results.iteration.variable(2).value(results.iteration.current);

disp(sprintf('data for time %f d at day %d', mod(tdata, 3), round(tdata/24)));

figure(100);
print('-depsc', [basedir '\results\temp\figure7-matlab.eps']);

