function dotheTestPlots()
% Comparison time plot, how do things look for KD 1nM vs KD 4uM

load resultsWithTolerableDrug.mat
plotit(AT_M_vector,t_record)
clear

hold on
load resultsWithInTolerableDrug.mat
plotit(AT_M_vector,t_record)
legend('K_D=4000nM','K_D=1nM','Location','SouthWest')
ylabel('#Proliferated Cells (DC+NaiveT+ActT)')
% % % % % 

clear
figure
load resultsWithTolerableDrugjustactT.mat
plotit(AT_M_vector,t_record)
clear

hold on
load resultsWithInTolerableDrugjustactT.mat
plotit(AT_M_vector,t_record)
legend('K_D=4000nM','K_D=1nM','Location','SouthWest')
ylabel('#Proliferated Cells (ActT only)')

function plotit(AT_M_vector,t_record)
% figure
LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
xlabeltext = 'Time (days)';
ylabeltext = '#Proliferated cells';
% legenddata = {'T-Cells (Epitope 1)','T-Cells (Epitope 2)'};

plot(t_record,AT_M_vector(:,1),'LineWidth',LT);

set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
% legend(legenddata,'Location','SouthEast');
xlabel(xlabeltext);
ylabel(ylabeltext);
set(gca,'yscale','log')
axis square

