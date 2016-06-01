function hColorBar = drawNormalizedVelocity (haxes, results, max_min, scale)

% draw the plot for the results

% set the figure axes to the surrenct axes
axes(haxes);

% get the system size (TODO later)
tubeRadius = 3500;


% padd with borders?
paddwithborders = false;

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = sprintf('%08d', itern); %string

% get either the solute or solid concentration
fn = results.directory;

variable = results.flow.available{results.flow.current};
fn = [results.directory '/flow/iteration' iteration '/' variable '.txt'];


if (nargin > 2)
    data = load(fn, 'ascii')*scale;
else
    data = load(fn, 'ascii');
end;


%get also the biomass
biomass = [];
for i = 1:length(results.solids.available)
    variable = results.solids.available{i};
    fn = [results.directory '/solids/iteration' iteration '/' variable '.txt'];
    temp = load(fn, 'ascii');
    if isempty(biomass),
        biomass = temp;
    else
        biomass = biomass + temp;
    end;
end;




gridSide = size(data, 1);
sys.X = tubeRadius / (0.5 - 1 / gridSide) * 1.01;
sys.Y = sys.X;

% dont forget to remove the extra 10 matrix entries from the level set matrix
if not((results.solids.show) |...
        (results.solutes.show) |...
        (results.flow.show)),
    data = data(10:end, :);
end;

%padd with boundary layer borders:
if (paddwithborders)
    pdata = paddWithBorders(data);
    pbiomass = paddWithBorders(biomass);
else
    %tube reactor case
    pdata = paddWithZeroFluxBorders(data);
    pbiomass = paddWithZeroFluxBorders(biomass);
end
[m, n] = size(data);
xgrid = ((0 : n+1) - 0.5) / (n) * sys.X;
ygrid = ((m+1 :-1: 0) - 0.5) / (m) * sys.Y;

[xgrid, ygrid] = meshgrid(xgrid, ygrid);
pbiomass((xgrid-sys.X/2).^2 + (ygrid-sys.Y/2).^2 >= tubeRadius.^2) = 10;
%pdata(pbiomass > 0) = -1;

% vector of contour lines to plot
minval = 0;
maxval = max(data(:));

%
lines = minval : (maxval-minval)/50 : maxval;

% draw the velocity field
[c,h] = contourf(xgrid, ygrid, pdata, lines);
colormap(jet(length(lines)));
%remove this line to get the black lines back
for i = 1:length(h),
    set(h(i), 'LineStyle', 'none');
end;

if (~results.particles),
    % draw the biomass field
    hold on;
    [C, h] = contour(xgrid, ygrid, pbiomass, 1e-4);
    hold off;
    %remove this line to get the black lines back
    for i = 1:length(h),
        set(h(i), 'LineColor', [1 1 1], 'LineWidth', 2);
    end;
end;


% set the color limits if requested
if (nargin > 2)
    if (~isempty(max_min))
        set(gca, 'CLimMode', 'manual', 'CLim', max_min);
    end;
end;

% pcolor(xgrid, ygrid, pdata);

axis equal;
axis xy; 
%set(gca, 'YLim', [-sys.Y*0.01 sys.Y], 'YTick', [0 :sys.Y*.1: sys.Y*.4, sys.Y]);
set(gca, 'YLim', [-sys.Y*0.01 sys.Y], 'YTick', []);
set(gca, 'XLim', [0 sys.X], 'XTick', []);
set(gca, 'FontSize', 8);


% colormap bone;
posAxes = get(gca, 'Position');
pcb(1) = (posAxes(1) + posAxes(3)) * 0.85;
pcb(2) = posAxes(2);
pcb(3) = (1-pcb(1))/5;
pcb(4) = posAxes(4);
% hcb = axes('position',pcb);
% colorbar(hcb);
hColorBar = colorbar();
%set(hcb, 'Position', pcb);

% set the figure axes to the surrenct axes again
axes(haxes);

% draw a Circle of the tube
center = [sys.X*0.5 sys.X*0.5];


% draw the particles if selected
if (results.particles),
    fn = [results.directory '/particles/iteration' iteration '.txt'];
    bacteria = load(fn, 'ascii');
    
    %draw the particles
    [m, n] = size(bacteria);    
    for i = 1:m,
        x = bacteria(i, 2);
        y = bacteria(i, 1);
        rcore = bacteria(i, 4);
        color = [bacteria(i, 5)/255, bacteria(i, 6)/255, bacteria(i, 7)/255]; 
        %color = [0.9 0.1 0.1]; 
        %color = [1 1 1]; 
        rcapsule = bacteria(i, 8);
        colorcapsule =...
            [bacteria(i, 9)/255, bacteria(i, 10)/255, bacteria(i, 11)/255]; 
        %draw capsule
        rectangle('Curvature', [1 1], 'Position',...
            [x-rcapsule y-rcapsule 2*rcapsule 2*rcapsule],...
            'FaceColor', colorcapsule, 'LineStyle', 'none');
        %draw core
        if (rcore > 0),
            rectangle('Curvature', [1 1], 'Position',...
                [x-rcore y-rcore 2*rcore 2*rcore],...
                'FaceColor', color, 'LineStyle', 'none');
        end;
    end;
end;

rectangle('Curvature', [1 1], ...
    'Position',...
    [center(1)-tubeRadius center(2)-tubeRadius 2*tubeRadius 2*tubeRadius], ...
    'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 3,...
    'EdgeColor', [1 1 1]);
