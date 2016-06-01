function growthCurveData = readGrowthCurve(fileName)

% read the growth curve datat into a data structure

if (isempty(dir(fileName)))
    warning(['file ' fileName ' not found']);
    return;
end;

% parse the numeric data all at once
fileData = importdata(fileName);

growthCurveData.time = fileData.data(:, 2);
growthCurveData.median = fileData.data(:, 21);
growthCurveData.average = fileData.data(:, 22);
growthCurveData.std = fileData.data(:, 23);
for i = 1:9,
    growthCurveData.sensor(i).value = fileData.data(:, 2 + i);
    growthCurveData.sensor(i).error = fileData.data(:, 11 + i);
end;
