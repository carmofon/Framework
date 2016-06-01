function drawConcentrationCurves (haxes, results)

% draw the plot for the results

% set the figure axes to the surrenct axes
axes(haxes);

% get the system size (TODO later)
sys.X = 500;
sys.Y = sys.X;

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = num2str(itern); %string

% get either the solute or solid concentration
fn = results.directory;

if (results.solids.show)
    variable = results.solids.available{results.solids.current};
    fn = [fn '/solids/iteration' iteration '/' variable '.txt'];
elseif (results.solutes.show)
    variable = results.solutes.available{results.solutes.current};
    fn = [fn '/solutes/iteration' iteration '/' variable '.txt'];
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

data = load(fn, 'ascii');

% dont forget to remove the extra 10 matrix entries from the level set matrix
if not((results.solids.show) | (results.solutes.show)),
    data = data(10:end, :);
end;

%padd with boundary layer borders:
pdata = paddWithBorders(data);
[m, n] = size(data);
 xgrid = ((0 : n+1) - 0.5) / (n) * sys.X;
 ygrid = ((m+1 :-1: 0) - 0.5) / (m) * sys.Y;
%xgrid = ((0:n+1) - 1) / (n) * sys.X;
%ygrid = ((m+1 :-1: 0) ) / (m) * sys.Y;


% vector of contour lines to plot
minval = 0;
maxval = max(data(:));

%
if ((results.solids.show) | (results.solutes.show))
    lines = minval : (maxval-minval)/10 : maxval;
else
    endint = tstep*100;
    step = log(endint-minval)/10;
    lines = [minval,  exp(step : step : log(endint))];
%     lines = [0,tstep];
end;

if isempty(lines),
    pcolor(xgrid, ygrid, pdata);
else,
    hold on;
    [c,h] = contour(xgrid, ygrid, pdata, lines, 'k-');
    hlab = clabel(c,h);
    hold off;
end;
