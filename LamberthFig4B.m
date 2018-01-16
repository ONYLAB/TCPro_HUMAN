seqtable = readtable('runpeptides.dat');

for i = 1:height(seqtable)
    
    seqid = seqtable{i,'SEQNAME'}{1};
    epitope = seqtable{i,'SEQUENCE'}{1};
    disp(seqid);
    %     isHuman(i) = isAntigenHuman(epitope);
    if mod(i,2)==0 %isHuman(i)==0 %Exogenous
        mkdir(seqid)
        cd(seqid)
        [responsesummary,significantresponsesummary,kon,responsevector,significancevector] = runalldonor(epitope);
        kons{i}=kon;
        allmeanresponses{i} = responsevector;
        allsignificance{i} = significancevector;
        responsesCPM(i)=responsesummary;
        responsesELI(i)=significantresponsesummary;
        cd ..
        save
    end
end

for i = 1:height(seqtable)
    if mod(i,2)==0 %isHuman(i)==0 %Exogenous
        for j = 1:51
            
            responses = allmeanresponses{1,i}{1,j}>=2;
            significance = allsignificance{1,i}{1,j}<=0.05;
            
            sigresponse(i,j,:) = responses.*significance;
            CUMsigresponses(i,j) = sign(sum(sigresponse(i,j,1:4)));
            
        end
    end
end

save