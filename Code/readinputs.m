function  [antigenname,cohortname,Nrun,MeanFp,StdFp,SampleConcentration,cohort,SIcutoff,sigcutoff] = readinputs()

%Read options file
fileID = fopen('options.dat');
options = textscan(fileID,'%s %*[^\n]');
fclose(fileID);

antigenname = options{1,1}{1,1};
cohortname = options{1,1}{2,1};
cohorthlafilename = options{1,1}{3,1};
fastaseq = options{1,1}{4,1};
Nrun = str2num(options{1,1}{5,1});
MeanFp = str2num(options{1,1}{6,1})/1e6;
StdFp = str2num(options{1,1}{7,1})/1e6;
SIcutoff = str2num(options{1,1}{8,1});
sigcutoff = str2num(options{1,1}{9,1});

locBLASTbin = options{1,1}{10,1};
locEpitopeDB = options{1,1}{11,1};
EpitopeDB  = options{1,1}{12,1};

SampleConcentration = read15mers(fastaseq,locBLASTbin,locEpitopeDB,EpitopeDB);

% %Temporary KDread
% % load KDeff.mat;
% load(['/home/osman.yogurtcu/Documents/NetMHCIIPanVersions/netMHCIIpan-3.2/Sequences/seqs/' antigenname '/Affinity.mat']);
% data = 1./sum(1./data);
% data = data(:);

Allabout15mers = readtable('15mers_KDs_Pepitopes.csv');
KDVals = Allabout15mers{:,2:end-1};
Pepitope = Allabout15mers{:,end};
data = 1./sum(1./KDVals,1); %Effective KD's per allele

%Read cohort HLAs
cohort = readtable(cohorthlafilename,'Delimiter',',');

%Find unique alleles
uniqHLA = unique([cohort{:,2};cohort{:,3}]);

%Assign allele indexes to the table    
t = readtable('allWorldAlleles.dat','Delimiter',',');
t.Donors=categorical(t.Donors);
for i = 1:length(uniqHLA)    
    
    cohort(strcmp(cohort{:,2},uniqHLA{i}),4) = {i};
    cohort(strcmp(cohort{:,3},uniqHLA{i}),5) = {i};
    KD = [];
    
    thehla = uniqHLA{i};
    thehla = ['DRB1_' thehla([end-4 end-3 end-1 end])];

    indexwithinAllWorldAlleles = find(t.Donors==thehla);
    KDeff = data(1,indexwithinAllWorldAlleles);
    
    cohort(strcmp(cohort{:,2},uniqHLA{i}),6) = {KDeff};
    cohort(strcmp(cohort{:,3},uniqHLA{i}),7) = {KDeff};
    
end
cohort(find(cohort{:,6}==cohort{:,7}),8)={1}; %If homozygous = 1

for donindex = 1:height(cohort)
  
    hla1 = strrep(strrep(cohort{donindex,2}{:}(5:end),'*','_'),':','');
    hla2 = strrep(strrep(cohort{donindex,3}{:}(5:end),'*','_'),':','');

    a1 = find(strcmp(hla1,string(t{:,1})));
    a2 = find(strcmp(hla2,string(t{:,1})));
    cohort(donindex,9) = {dot(sum(1./KDVals(:,[a1 a2]),2),Pepitope)};
    cohort(donindex,10) = {cohort{donindex,9}/sum([1/cohort{donindex,6} 1/cohort{donindex,7}])};
    
end

cohort.Properties.VariableNames(4) = {'Allele1index'};
cohort.Properties.VariableNames(5) = {'Allele2index'};
cohort.Properties.VariableNames(6) = {'KDeff1'};
cohort.Properties.VariableNames(7) = {'KDeff2'};
cohort.Properties.VariableNames(8) = {'Homozygous'};
cohort.Properties.VariableNames(9) = {'PhiNumerator'};
cohort.Properties.VariableNames(10) = {'Phi'};

movefile('15mers.csv','../Output/');
movefile('15mers_KDs_Pepitopes.csv','../Output/');
movefile('15mers_Ranks.csv','../Output/');
movefile('NetMHCIIpan_out.csv','../Output/');
delete('input.fa');
delete('results.txt');
delete('out.dat');
delete('NetMHCIIpan_out.xls')

cd ../Output
% Write cohort member Effective KDs into output folder
writetable(cohort,[antigenname '_' cohortname '_cohortKDeff.csv'],'Delimiter',',');
cd ../Input