function roughness = plotRoughness(dirName)

% results = getResultsFromDirectory(dirName)
% Parse directory structure for results

results = getResultsFromDirectory(dirName);

% get the names of solutes from the first iteration
iterationList = dir([dirName '/solids/iteration*']);
roughness = zeros(1, length(iterationList));
for i = 1:length(iterationList),
    itern = i-1; %number
    iteration = sprintf('%08d', itern); %string
    variable = results.solids.available{results.solids.current};
    fn = [results.directory '/solids/iteration' iteration '/' variable '.txt'];
    data = load(fn, 'ascii');
    if i ==1
        [m,n] = size(data);
        m = m-(1:m);
        n = 1:n;
        [n, m] = meshgrid(n,m);
    end
    data(data>0) = m(data>0);
    roughness(i) = var(max(data, [], 1));
end;

plot(roughness);
