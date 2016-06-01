% scMovieMaker2

basedir = 'D:\results\';

data(1).dirName = [basedir 'drg-detachment-eps\umax5.47E-01kdet1.00E-03kd_eps1.40E-03ki1.00E-03'];
data(2).dirName = [basedir 'drg-detachment-eps-limitation\umax5.47E-01kdet1.00E-03keps5.00E01ki1.00E-03'];
data(3).dirName = [basedir 'drg-detachment-AI-eps\umax5.47E-01kdet1.00E-03rai5.00E-02kai3.50E-04'];

% for i = 1:length(data),
i = 3;
makeMovie(data(i).dirName);
% end;




