function plotresponse(Day,TYPE,n)
% Time dependent Cell Count plot

if TYPE == 'S'
    load(['D' num2str(Day) '_SAMPLE_n' num2str(n) 'results.mat']);
else
    load(['D' num2str(Day) '_CULTURE_n' num2str(n) 'results.mat']);
end

% % % % % % % % % % % % % % % % %
figure
LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
xlabeltext = 'Time (days)';
ylabeltext = '#cells';
legenddata = {'Immature DC','Mature DC','Naive T','Activated T','Bystander T'};

plot(t_record,ID_vector(:,1),'LineWidth',LT);
hold on
plot(t_record,MD_vector(:,1),'LineWidth',LT);
plot(t_record,NT_vector(:,1),'LineWidth',LT);
plot(t_record,AT_N_vector(:,1),'LineWidth',LT);
plot(t_record,MT_vector(:,1),'LineWidth',LT);

set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
legend(legenddata,'Location','eastoutside');
xlabel(xlabeltext);
ylabel(ylabeltext);
set(gca,'yscale','log')
axis square