% scMovieMAker

basedir = 'D:\results\detachment\constantconcentration\dfreps\';

data(1).dirName = [basedir 'lixo'];
% data(1).dirName = [basedir 'Cs1e-3Kd2e-4'];
% data(2).dirName = [basedir 'Cs2e-3Kd1e-3'];
% data(3).dirName = [basedir 'Cs2e-3Kd5e-3'];
% data(4).dirName = [basedir 'Cs4e-3Kd1e-2'];
% data(5).dirName = [basedir 'Cs4e-3Kd5e-3'];

for i = 1:length(data),
    makeMovie(data(i).dirName);
end;




