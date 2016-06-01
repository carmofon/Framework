clear
%close all

cd /Users/CarmonaFontaine/Documents/Xlab/Models/160503_Angiogenesis/
FolderUni = 'TestUni';
FolderTargeted = 'TestTargeted';

cd (FolderUni)
uni = readtable('simulationResults.txt','Delimiter','\t');
cd ..
cd (FolderTargeted)
targeted = readtable('simulationResults.txt','Delimiter','\t');
cd ..

%%
subplot(2,2,1)
plot((uni.activeC_M_))
hold on
plot((targeted.activeC_M_))
hold off
legend('Uniform','Targeted','location','northwest')
set(gca,'yscale','lin')
title ('cell growth')
xlabel('Time')

subplot(2,2,2)
plot((uni.VEGF_M_))
hold on
plot((targeted.VEGF_M_))
hold off
set(gca,'yscale','lin')
title ('total VEGF')
xlabel('Time')