function data = readCurrentSpeciesData(results)

% get the current iteration number
itern = results.iteration.variable(1).value(results.iteration.current); %number
iteration = num2str(itern);                                             %string

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
