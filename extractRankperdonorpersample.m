function extractRankperdonorpersample()
% This function plots rank or KD info per donor per drug

close all
seqtable = readtable('runpeptides.dat');
colid = 'bbggrrcc';
list = [1 5 2 6 3 7 4 8];

for i = 1:height(seqtable)
    
    seqid = seqtable{i,'SEQNAME'}{1};
    cd(seqid)
    subplot(2,4,list(i));
    for d = 1:51
       cd(num2str(d))
        netmhcIItbl = readtable('out.dat');
        rank(i,d,:) = netmhcIItbl{1,[6 9]};
        kd(i,d,:) = netmhcIItbl{1,[6 9]-1};
        donrank(d,i) = mean(mean(rank(i,d,:)));
       cd ..
       plothla(d,rank(i,d,:),colid(i));
    end
    cd ..
%     ylabel('K_D (nM)');
    ylabel('%Affinity rank');
    xlabel('DonorID');
    set(gcf,'color','w');
    title(seqid)
    axis square
    xlim([0 52])
    avgrank(i) = mean(mean(rank(i,:,:)));
    save('ranks.mat');
    
%     h=histogram(donrank(:,i),round(sqrt(51)),'Normalization','probability');    
%     edges = h.BinEdges;
%     width = h.BinWidth;
%     centers(i,:) = edges(1:end-1) + width/2;
%     vals(i,:) = h.Values;
%     close
    
end

% figure
% axis square
% set(gcf,'color','w');
% xlabel('K_D (nM)');
% ylabel('%Rank Affinity')

% LT = 3; %Line thickness
% AxFS = 24; %Ax Fontsize
% AxLW = 2; %Ax LineWidth
% xlabeltext = '%Average pMHC Affinity Rank';
% ylabeltext = 'Probability over 51 donors';
% plot(centers',vals','LineWidth',LT)
% set(gcf,'color','w');
% set(gca,'fontsize', AxFS);
% set(gca,'LineWidth',AxLW);
% xlabel(xlabeltext);
% ylabel(ylabeltext);
% axis square
% legend(seqtable{:,'SEQNAME'},'Location','NorthWest')

function plothla(donorind,data,colid)

plot([donorind donorind],data(:),colid)
hold on
scatter([donorind donorind],data(:),colid,'.')