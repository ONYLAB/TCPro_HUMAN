function [responsevector,decision] = runalldonor(theepitope,ParChangeIndex,ParChange)

allelelist = readtable('detailedAlleleshaplotypeddonorNETMHCIIreadablewithoutDQBno345.dat','Delimiter',',');
epitopes{1} = theepitope;
mypwd = pwd;
addpath(mypwd);
tic
for donor_ID = 1:1%height(allelelist)
    
    disp(['Donor: ' num2str(donor_ID)]);
    
    [s,m] = unix(['rm -rf ' num2str(donor_ID)]);
    mkdir(num2str(donor_ID));
    cd(num2str(donor_ID));
    HLA_DR{1} = allelelist{donor_ID,1}{1};
    HLA_DR{2} = allelelist{donor_ID,2}{1};
    [responsevector(donor_ID,:,:),decision(1,donor_ID)] = run3samplesDaylimit(epitopes,HLA_DR,donor_ID,ParChangeIndex,ParChange);
    cd ..
%     save
end
toc
% cpmresponse = 100*sum(maxcpmresponsevector>2)/height(allelelist);
% eliresponse = 100*sum(maxeliresponsevector>2)/height(allelelist);

function [responsevector,decision] = run3samplesDaylimit(epitopes,HLA_DR,donor_ID,ParChangeIndex,ParChange)

colnames = {};
colnameindex = 0;
coldata = [];
n=3;
for Daylimit = 5:8
    for i = 1:n
        [response(Daylimit-4,i,1),kon,ELISPOT(Daylimit-4,i,1)] = Main_human(Daylimit,0,epitopes,HLA_DR,donor_ID*i,ParChangeIndex,ParChange);%SimType=0, with sample
        status = movefile('Parameters.mat',['D' num2str(Daylimit) '_CULTURE_n' num2str(i) 'Parameters.mat']);
        status = movefile('results.mat',['D' num2str(Daylimit) '_CULTURE_n' num2str(i) 'results.mat']);
        
        [response(Daylimit-4,i,2),kon,ELISPOT(Daylimit-4,i,2)] = Main_human(Daylimit,1,epitopes,HLA_DR,donor_ID*i,ParChangeIndex,ParChange);%SimType=1, with sample
        status = movefile('Parameters.mat',['D' num2str(Daylimit) '_SAMPLE_n' num2str(i) 'Parameters.mat']);
        status = movefile('results.mat',['D' num2str(Daylimit) '_SAMPLE_n' num2str(i) 'results.mat']);
    end

    responsevector(Daylimit-4,:) = response(Daylimit-4,:,2)./response(Daylimit-4,:,1);
end

decision = sign(sum(sum(responsevector>1.9)));

% ELISPOTresp = mean(ELISPOT(Daylimit-4,:,2)) / mean(ELISPOT(Daylimit-4,:,1));
% [~,p,~,~] = ttest2(ELISPOT(Daylimit-4,:,2),ELISPOT(Daylimit-4,:,1));
% significancevector = [significancevector; p];
% responsevector = [responsevector; ELISPOTresp];
% 
% responsesummary = sum(responsevector>2); %out of 5
% significantresponsesummary = sum((responsevector>2).*(significancevector<0.05));

% sem = stdresp;
% LT = 3; %Line thickness
% AxFS = 24; %Ax Fontsize
% AxLW = 2; %Ax LineWidth
% xlabeltext = 'Time (days)';
% ylabeltext = 'Response: Stimulation Index';
% errorbar(5:8,meanresp,sem,'LineWidth',LT);
% set(gcf,'color','w');
% set(gca,'fontsize', AxFS);
% set(gca,'LineWidth',AxLW);
% xlabel(xlabeltext);
% ylabel(ylabeltext);
% axis square
% title(['Donor#' num2str(donor_ID)]);
% save([num2str(donor_ID) '.mat'],'response','ELISPOTresp','responsevector','significancevector');
% close

% T = table(coldata','RowNames',colnames);
% writetable(T,'myDonorData.dat','WriteRowNames',true);

% maxmeanresp = max(meanresp);
% maxmeanresp = max(theresponse(:,2))/max(theresponse(:,1));