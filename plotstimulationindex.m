load

LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
xlabeltext = 'Donor ID';
ylabeltext = 'Stimulation Index';

for i =1:51
    
   errorbar(i,meanresp{1,i}(4),stdresp{1,i}(4),'LineWidth',LT);
   hold on
    
end

set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
% legend(legenddata,'Location','SouthEast');
xlabel(xlabeltext);
ylabel(ylabeltext);

axis square
title('158V')
xlim([0 52])
ylim([0 5])
plot(0:52,2*ones(length(0:52),1),'k:','LineWidth',3)