function donres = plotscatterdonor()
% This function is for showing the stimulation index values per donor per
% drug
close all
load('matlab.mat','allmeanresponses','seqtable');

for i = 1:height(seqtable)
    
    seqid = seqtable{i,'SEQNAME'}{1};
    epitope = seqtable{i,'SEQUENCE'}{1};
    
    responsevector = allmeanresponses{i};
    
    figure
    for j = 1:51
        
        donordata = responsevector{j};
        plotindivscatter(j,donordata);
        donres(i,j) = sum(donordata(1:4)>2);
    end
    title(seqid)
    plot(1:51,2*ones(51,1),'--','LineWidth',3)
    ylim([0 5])
    xlim([0 52])
end

function plotindivscatter(donorid,data)

c = linspace(1,4,4);
sz = 25;
% for i = 1:4
scatter(donorid*ones(4,1),data(1:4),sz,c,'filled');
hold on
% end
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
% axis square
ylabel('Stimulation Index')
xlabel('Donor Index')
