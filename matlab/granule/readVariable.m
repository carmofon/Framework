function data2D = readVariable(results)

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = sprintf('%08d', itern); %string

% get either the solute or solid concentration
fn = results.directory;

if (results.solids.show);
    variable = results.solids.available{results.solids.current};
    fn = [fn '/solids/iteration' iteration '/' variable '.txt'];
else (results.solutes.show);
    variable = results.solutes.available{results.solutes.current};
    fn = [fn '/solutes/iteration' iteration '/' variable '.txt'];
end;

%
data2D = load(fn, 'ascii');
