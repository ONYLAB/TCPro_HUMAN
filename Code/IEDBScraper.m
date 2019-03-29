function KD = IEDBScraper(peptidesequence,nallele,allelelist)

url = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/';
method = 'method';
methodSelect = 'NetMHCIIpan';
sequence_text = 'sequence_text';
sequence_textSelect = peptidesequence;
allele = 'allele';
alleleSelect = '';
for i = 1:nallele
    if i>1
        alleleSelect = [alleleSelect ','];
    end
    al = allelelist{i};
    alleleSelect = [alleleSelect al]; % alleleSelect = 'HLA-DRB1*10:01,HLA-DRB1*12:01';
end

count = 0;
err_count = 0;
while count == err_count
    try
        response = webwrite(url,method,methodSelect,sequence_text,sequence_textSelect,allele,alleleSelect,weboptions('Timeout',24*60*60));
    catch MyErr
        err_count = err_count + 1;
        disp('IEDB Connection failed, retrying in 1 minute');
        pause(60);
    end
    count = count + 1;
end

%Scape KD data from IEDB per 15-mers
lines = splitlines(response);
for i = 2:length(lines)-1
    temp = sscanf(lines{i},'%s%i%i%i%s%s%f%f');
    KD(i-1,1) = temp(end-1);
end