function [responsevector,decision] = runalldonor(theepitope,randindex)

allelelist = readtable('detailedAlleleshaplotypeddonorNETMHCIIreadablewithoutDQBno345.dat','Delimiter',',');
epitopes{1} = theepitope;
mypwd = pwd;
addpath(mypwd);
tic
for donor_ID = 1:height(allelelist)
    
    disp(['Donor: ' num2str(donor_ID)]);
    
    [s,m] = unix(['rm -rf ' num2str(donor_ID)]);
    mkdir(num2str(donor_ID));
    cd(num2str(donor_ID));
    delete('*.mat')
    HLA_DR{1} = allelelist{donor_ID,1}{1};
    HLA_DR{2} = allelelist{donor_ID,2}{1};
    [responsevector(donor_ID,:,:),decision(1,donor_ID)] = run3samplesDaylimit(epitopes,HLA_DR,donor_ID*randindex);
    cd ..
    %     save
end
toc
% cpmresponse = 100*sum(maxcpmresponsevector>2)/height(allelelist);
% eliresponse = 100*sum(maxeliresponsevector>2)/height(allelelist);

function [responsevector,decision] = run3samplesDaylimit(epitopes,HLA_DR,donor_ID)

n=3;
Daylimit = 8;
for i = 1:n
    [response0(i,:),kon,ELISPOT(i,1)] = Main_human(Daylimit,0,epitopes,HLA_DR,donor_ID*i);%SimType=0, with sample
    status = movefile('Parameters.mat',['D' num2str(Daylimit) '_CULTURE_n' num2str(i) 'Parameters.mat']);
    status = movefile('results.mat',['D' num2str(Daylimit) '_CULTURE_n' num2str(i) 'results.mat']);
    
    [response1(i,:),kon,ELISPOT(i,2)] = Main_human(Daylimit,1,epitopes,HLA_DR,donor_ID*i);%SimType=1, with sample
    status = movefile('Parameters.mat',['D' num2str(Daylimit) '_SAMPLE_n' num2str(i) 'Parameters.mat']);
    status = movefile('results.mat',['D' num2str(Daylimit) '_SAMPLE_n' num2str(i) 'results.mat']);
end

responsevector = response1./response0;
for days = 1:4
    
    [~,sigvector(1,days),~,~] = ttest2(response1(:,days),response0(:,days));
    
end

divisioncriterion = 1.9;
sigcriterion = 0.05;
decisionperday = sum(responsevector>divisioncriterion).*(sigvector<sigcriterion);
decision = sign(sum(decisionperday));
