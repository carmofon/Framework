function hColorBar = drawResults (haxes, results, max_min, scale)

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

if (results.solids.show)
    variable = results.solids.available{results.solids.current};
    fn = [fn '/solids/iteration' iteration '/' variable '.txt'];
elseif (results.solutes.show)
    variable = results.solutes.available{results.solutes.current};
    fn = [fn '/solutes/iteration' iteration '/' variable '.txt'];
elseif (results.flow.show)
    variable = results.flow.available{results.flow.current};
    fn = [fn '/flow/iteration' iteration '/' variable '.txt'];
else 
    % for showing the density level set
    fn = [fn '/detachmentLevelSet/iteration' iteration '.txt'];
    % get the value of the time step
    if (itern > 1)
        tstep = results.iteration.variable(2).value(itern)...
            - results.iteration.variable(2).value(itern-1);
    else
        tstep = results.iteration.variable(2).value(itern);
    end;
end;

data = load(fn, 'ascii')*scale;
%
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
else
    pdata = paddWithZeroFluxBorders(data);
end
[m, n] = size(data);
xgrid = ((0 : n+1) - 0.5) / (n) * sys.X;
ygrid = ((m+1 :-1: 0) - 0.5) / (m) * sys.Y;


% vector of contour lines to plot
minval = 0;
maxval = max(data(:));

%
if ((results.solids.show) | (results.solutes.show) | (results.flow.show))
    lines = minval : (maxval-minval)/50 : maxval;
else
    endint = tstep*100;
    step = log(endint-minval)/10;
    lines = [minval,  exp(step : step : log(endint))];
%     lines = [0,tstep];
end;

if isempty(lines),
    pcolor(xgrid, ygrid, pdata);
else,
    [c,h] = contourf(xgrid, ygrid, pdata, lines);
    colormap(jet(length(lines)));
%     newmap = bone(round(length(lines)*2)) ;
%     newmap = newmap(end-length(lines):end,:);
%     colormap(newmap(end:-1:1,:));
    %remove this line to get the black lines back
    for i = 1:length(h), 
        set(h(i), 'LineStyle', 'none');
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

rectangle('Curvature', [1 1], ...
    'Position',...
    [center(1)-tubeRadius center(2)-tubeRadius 2*tubeRadius 2*tubeRadius], ...
    'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 3,...
    'EdgeColor', [0.9 0.1 0.1]);


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
        %color = [bacteria(i, 5)/255, bacteria(i, 6)/255, bacteria(i, 7)/255]; 
        %color = [0.9 0.1 0.1]; 
        color = [1 1 1]; 
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
    
    % draw the substratum circle
    %     rectangle('Position', [0 -sys.Y*0.01 sys.Y sys.Y*0.01],...
    %         'FaceColor', [0.3, 0.3, 0.3]);
end;
