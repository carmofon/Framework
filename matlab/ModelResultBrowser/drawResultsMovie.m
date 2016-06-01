function drawResultsMovie (haxes, results)

% draw the plot for the results

% set the color to white 
set(gcf, 'Color', [1 1 1]);

% set the figure axes to the surrenct axes
axes(haxes);

% get the system size (TODO later)
sys.X = 2000;
sys.Y = sys.X;

% get either the solute or solid concentration
fn = results.directory;
if (results.solids.show)
    variable = results.solids.available{results.solids.current};
    fn = [fn '/solids/iteration'...
            num2str(results.iteration.available(results.iteration.current))...
            '/' variable '.txt'];
else 
    variable = results.solutes.available{results.solutes.current};
    fn = [fn '/solutes/iteration'...
            num2str(results.iteration.available(results.iteration.current))...
            '/' variable '.txt'];
end;

data = load(fn, 'ascii');

%padd with boundary layer borders:
pdata = paddWithBorders(data);
[m, n] = size(data);
% xgrid = ((0 : n+1) - 0.5) / (n) * sys.X;
% ygrid = ((m+1 :-1: 0) - 0.5) / (m) * sys.Y;
xgrid = ((0:n+1) - 1) / (n) * sys.X;
ygrid = ((m+1 :-1: 0) ) / (m) * sys.Y;

% vector of contour lines to plot (linear)
% minval = min(data(:));
% maxval = max(data(:));
% lines = minval : (maxval-minval)/10 : maxval;
% vector of contour lines to plot (log)
%minval = exp(min(data(:)));
minval = exp(0);
maxval = exp(max(data(:)));
lines = log(minval : (maxval-minval)/10 : maxval);

[c,h] = contourf(xgrid, ygrid, pdata, lines);
%pcolor(xgrid, ygrid, pdata);

newmap = bone(round(length(lines)*2)) ;
newmap = newmap(end-length(lines):end,:);
colormap(newmap(end:-1:1,:));
colorbar;

% draw the particles if selected
if (results.particles),
    fn = [results.directory '/particles/iteration'...
            num2str(results.iteration.available(results.iteration.current))...
            '.txt'];
    bacteria = load(fn, 'ascii');
    
    %draw the particles
    [m, n] = size(bacteria);    
    for i = 1:m,
        x = bacteria(i, 4);
        y = bacteria(i, 3);
        r = bacteria(i, 1);
        rcore = bacteria(i, 2);
        color = [bacteria(i, 6)/255, bacteria(i, 7)/255, bacteria(i, 8)/255]; 
        rectangle('Curvature', [1 1], 'Position',...
            [x-r y-r 2*r 2*r], 'FaceColor', 'none', 'LineStyle', '-');
        rectangle('Curvature', [1 1], 'Position',...
            [x-rcore y-rcore 2*rcore 2*rcore], 'FaceColor', color, 'LineStyle', 'none');
    end;
    
    % draw the substratum rectangle
    rectangle('Position', [0 -sys.Y*0.01 sys.Y sys.Y*0.01], 'FaceColor', [0.3, 0.3, 0.3]);    
end;
axis equal;
axis xy; 
set(gca, 'YLim', [-sys.Y*0.01 sys.Y]);
set(gca, 'XLim', [0 sys.X]);
hold off;
