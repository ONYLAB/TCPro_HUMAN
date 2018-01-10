function isHuman = isAntigenHuman(seq) 

isHuman = 0;
Antigen.Sequence = seq; %eg 'SKPQGRIVGGKDCPKGECPWQVL'
[RID1,ROTE] = blastncbi(Antigen,'blastp','Database','swissprot','Entrez','Human[organism]','MaxNumberSequences',20);
report1 = getblast(RID1,'WaitTime',ROTE);

for i = 1:length(report1.Hits)
    HumanSequence = report1.Hits(i).Hsps(1).Alignment(3,:);
    if strcmp(HumanSequence,Antigen.Sequence)
        isHuman = 1;
        disp('Antigen sequence is endogenous');
        break;
    end
end