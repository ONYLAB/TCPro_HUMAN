function RankvsnmperAlleleperseqlength()
% This function is reading in NetMHCIIpan thresholds data in Logk50 units 
% and converting them into nM values

conversiontable = readtable('pseudosequences.dat');
mypwd = pwd;

alleles = conversiontable{:,1};

for i = 1:height(conversiontable)
    
    seq = conversiontable{i,1};
    if strcmp('DRB',seq{1}(1:3))
        %YMFFQFSGGAILNTLFGQFEYFDLEKVRVHLGMT.thr_cmb
        pseudoseq = [pwd '\thresholds\' conversiontable{i,2}{1} '.thr_cmb'];
        pseudoseqintxt = [pwd '\thresholds\' conversiontable{i,2}{1} '.txt'];
        status = movefile(pseudoseq,pseudoseqintxt);
        convdata = readtable(pseudoseqintxt,'HeaderLines',1);
        convdata = convdata(:,2:end);%Format: thr 9 10 11 12 13 14 15 16 17 18 19 (PepSeqLength)
        conversiondata(i,:,:)= 5e4.^(1-convdata{:,2:end}); %IC;
    end
end

DRBalleles = alleles(1:length(conversiondata),1);

save('ConversionDatarankvsnm.mat','DRBalleles','conversiondata');