seqtable = readtable('runpeptides.dat');

for i = 1:height(seqtable)
   
    seqid = seqtable{i,'SEQNAME'}{1};
    epitope = seqtable{i,'SEQUENCE'}{1};
    mkdir(seqid)
    cd(seqid)
    [cpmresponse,eliresponse,kon,responsevector,significancevector] = runalldonor(epitope);
    kons{i}=kon;
    allmeanresponses{i} = responsevector;
    allsignificance{i} = significancevector;
    responsesCPM(i)=cpmresponse;
    responsesELI(i)=eliresponse;
    cd ..
    save
end