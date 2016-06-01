function drawPartilcesInGrayscale (results)

% draw the plot for the results

%
sys.X = results.sizeX;
sys.Y = sys.X;

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = sprintf('%08d', itern); %string

% get either the solute or solid concentration
fn = results.directory;




% draw the particles if selected
fn = [results.directory '/particles/iteration' iteration '.txt'];
bacteria = load(fn, 'ascii');

%draw the particles
[m, n] = size(bacteria);
for i = 1:m,
    x = bacteria(i, 2);
    y = bacteria(i, 1);
    rcore = bacteria(i, 4);
    color = [bacteria(i, 5)/255, bacteria(i, 6)/255, bacteria(i, 7)/255];
    if or(isred(color), isblue(color))
        color = [0.2 0.2 0.2];
    elseif isgreen(color)
        color = [1 1 1];
    else
        color = [0.5 0.5 0.5];
    end;
    
    
    %draw core
    rectangle('Curvature', [1 1], 'Position',...
        [x-rcore y-rcore 2*rcore 2*rcore],...
        'FaceColor', color, 'LineStyle', 'none');
end;

axis equal tight;
axis xy;
set(gca, 'YLim', [0 sys.Y], 'YTick', []);
set(gca, 'XLim', [0 sys.X], 'XTick', []);
set(gca, 'Color', [0 0 0]);


function a = isred(c)

a = and((c(1) > c(2)), (c(1) > c(3)));

function a = isblue(c)

a = and((c(3) > c(1)), (c(3) > c(2)));

function a = isgreen(c)

a = and((c(2) > c(1)), (c(2) > c(3)));