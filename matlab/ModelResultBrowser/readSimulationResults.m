function results = readSimulationResults(filename, resultsin)

% resultsout = readSimulationResults(filename, resultsin) 
% read the data in a simulationResults.txt file
%
% Joao Xavier (j.xavier@tnw.tudelft.nl) - February 2004
% improved using dlmread August 2005

results = resultsin;

if (isempty(dir(filename)))
    warning(['file ' filename ' not found']);
    return;
end;

% if the file exists, open it
fid = fopen(filename);

% first line is the header
tline = fgetl(fid);
n = [strfind(tline, sprintf('\t')), length(tline)];
%get the name of each variable
for i = 1:length(n),
    if i == 1,
        ini = 1;
    else,
        ini = n(i-1)+1;
    end;
    if (i == length(n)),
        fin = n(i);
    else,
        fin = n(i)-1;
    end;
    results.iteration.variable(i).name = tline(ini:fin);
end;

% old code:
% % parse the numeric data in file one line at a time
% j = 1;
% while 1
%     tline = fgetl(fid);
%     if ~ischar(tline), break, end
%     data = eval(['[' tline ']']);
%     for i = 1:length(n),
%         results.iteration.variable(i).value(j) = data(i);
%     end;
%     j = j+1;
% end
% fclose(fid);

%improved code!!
fclose(fid);
% parse the numeric data all at once
data = dlmread(filename, '\t', 1, 0);
for i = 1:length(n),
    results.iteration.variable(i).value = data(:, i);
end;

