figure
% load
LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
xlabeltext = 'Time (days)';
ylabeltext = '#Proliferated cells';
% legenddata = {'T-Cells (Epitope 1)','T-Cells (Epitope 2)'};

plot(t_record,AT_M_vector(:,1),'LineWidth',LT);
hold on
plot(t_record,MT_vector(:,1),'LineWidth',LT);
hold on
plot(t_record,FT_vector(:,1),'LineWidth',LT);
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
% legend(legenddata,'Location','SouthEast');
xlabel(xlabeltext);
ylabel(ylabeltext);
set(gca,'yscale','log')
axis square
legend('DC+NT+AT','DC','NT')