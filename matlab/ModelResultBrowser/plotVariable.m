function [] = plotVariable(results, varn, multiplier, haxes)

axes(haxes);

iteration = results.iteration.current;
times = results.iteration.variable(2).value./24; 
vars = multiplier*(results.iteration.variable(varn).value); 
hpl = plot(times(1:iteration), vars(1:iteration), 'r-');
set(hpl, 'LineWidth', 1);
hxl = xlabel('Time [day]');
set(hxl, 'FontSize', 8);
hyl = ylabel(results.iteration.variable(varn).name);
set(hyl, 'FontSize', 8);
set(haxes, 'YLim', [min(vars) max(vars)]);
set(haxes, 'XLim', [min(times) max(times)]);
set(haxes, 'FontSize', 5);