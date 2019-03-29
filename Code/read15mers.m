function Concentration = read15mers(fastaseq,locBLASTbin,locEpitopeDB,EpitopeDB)

ind = 0;
totalPlength = 0;
myseq = fastaread(fastaseq);
numchain = length(myseq);

x = 15; %15mers only
MaximumNumberofxmers = 0; %Max. num possible x-mers cleaved-out

for ns = 1:numchain
    
    chain = myseq(ns).Sequence;
    
    chains{ns,1} = chain;    
    chainID(ns) = ns;
    
    for k = 1:ns-1
        if strcmp(chains{ns,1},chains{k,1})
            chainID(ns) = chainID(k);
            break;
        end
    end
    
    L = length(chain);
    totalPlength = totalPlength + L;
    MaximumNumberofxmers = MaximumNumberofxmers + floor(L/x);
    for i = 8:length(chain)-7
        ind = ind+1;
        seqs{ind,1} = chain(i-7:i+7);
        seqs{ind,2} = ns;
        seqs{ind,3} = chainID(ns);
    end
    
end

if totalPlength >=39 %So that Exenatide is a protein-drug
    %classify as protein
    T = table((1:ind)',[seqs{:,2}]',[seqs{:,3}]',[seqs(:,1)]);
else
    %classify as peptide    
    T = table(1,1,1,{chain});    
end
T.Properties.VariableNames = {'Index','ChainNum','ChainID','xmers'};
writetable(T,'15mers.csv','Delimiter',',');

if totalPlength >=39 %So that Exenatide is a protein-drug
    %classify as protein
    PofCleavability = MaximumNumberofxmers/ind;
    Concentration = PofCleavability*0.3; %MicroMolar
else
    %classify as peptide
    Concentration = 5.0; %MicroMolar
end

%% Run NetMHCIIPan
HLAs = readtable('allWorldAlleles.dat','Delimiter',',');
HLAs = sortrows(HLAs);
HLAlist = '';
for i = 1:height(HLAs)
    HLAlist = [HLAlist ',' HLAs{i,1}{:}];
end
HLAlist = HLAlist(2:end);

NetMHCIIPanlocation = '/home/osman.yogurtcu/Documents/NetMHCIIPanVersions/netMHCIIpan-3.2/';

if totalPlength >=39
    % classify as protein
    command = [NetMHCIIPanlocation 'netMHCIIpan -inptype 0 -f sequence.fasta -a ' HLAlist ' -xls NetMHCIIpan_out.xls > out.dat'];
else
    % classify as peptide
    unix('awk ''{if (NR!=1) {print}}'' sequence.fasta > peptide.fasta');
    command = [NetMHCIIPanlocation 'netMHCIIpan -inptype 1 -f peptide.fasta -a ' HLAlist ' -xls NetMHCIIpan_out.xls > out.dat'];
end

if exist('../Output/NetMHCIIpan_out.csv') == 2 && totalPlength >=39  %TEMPORARY  
    status = unix('cp ../Output/NetMHCIIpan_out.csv .');    
else
    status = unix(command);
    status = unix('cp NetMHCIIpan_out.xls NetMHCIIpan_out.csv');
end
tbl = readtable('NetMHCIIpan_out.csv');
tbl(:,end-1:end) = [];
Kds = tbl(:,[2 5:3:end]);
Kds.Properties.VariableNames(2:end) = HLAs{:,1};

Ranks = tbl(:,[2 6:3:end]);
Ranks.Properties.VariableNames(2:end) = HLAs{:,1};

%% Obtain Probability of CD4+ Epitope Likelihood
if totalPlength >=39
    % classify as protein
    for i = 1:height(Kds)
        [foreign(i,1),~,~,~] = isAntigenHuman(Kds{i,1}{:},locBLASTbin,locEpitopeDB,EpitopeDB);
    end
else
    % classify as peptide
    foreign = 1.0;
end

Kds{:,width(Kds)+1} = foreign;
Kds.Properties.VariableNames(end) = {'P_epitope'};

writetable(Kds,'15mers_KDs_Pepitopes.csv','Delimiter',',');
writetable(Ranks,'15mers_Ranks.csv','Delimiter',',');


function [ProbForeign,Data,alignmentredundancymat,lenuniqscores] = isAntigenHuman(seq,locBLASTbin,locEpitopeDB,EpitopeDB)

myseq(1).Sequence = seq; %eg 'SKPQGRIVGGKDCPKGECPWQVL'
myseq(1).Header = '15mer';
if exist('input.fa')==2
    delete('input.fa');
end
fastawrite('input.fa',myseq);

k = 1;
a = 7.5;
ProbForeign = 0;
alignmentredundancymat = [];
lenuniqscores = nan;

command = [locBLASTbin 'blastp -ungapped -db ' locEpitopeDB EpitopeDB ' -query input.fa -out results.txt -sum_stats T -evalue 200000 -matrix PAM30 -comp_based_stats 0'];
status = unix(command);

Data = blastreadlocal('results.txt',0);

numhits = length(Data.Hits);

if numhits>0
    
    for i = 1:numhits
        
        positives = Data.Hits(i).HSPs(1).Positives.Match;
        identities = Data.Hits(i).HSPs(1).Identities.Match;
        scores(i) = identities + 0.5*(positives-identities);
        
    end
    
    alignmentredundancymat = zeros(numhits);
    for i = 1:numhits-1
        for j = i+1:numhits
            %disp([i j])
            subseqavail = length(strfind(Data.Hits(i).HSPs(1).Alignment(3,:),Data.Hits(j).HSPs(1).Alignment(3,:)));
%             fulseqavail = strcmp(Data.Hits(i).HSPs(1).Alignment,Data.Hits(j).HSPs(1).Alignment);
            alignmentredundancymat(i,j) = sign(subseqavail);
        end
    end
    
    scores(find(sum(alignmentredundancymat)>0))=[]; %Eliminate redundancy
    lenuniqscores = length(scores);
    
    ex = sum(exp(-k*(a-scores)));
    ProbForeign = ex/(1+ex);
    
end
