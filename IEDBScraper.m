function rank = IEDBScraper(peptidesequence,nallele,allelelist)

url = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/';
method = 'method';
methodSelect = 'recommended';
sequence_text = 'sequence_text';
sequence_textSelect = peptidesequence;
allele = 'allele';
alleleSelect = '';
for i = 1:nallele
    if i>1
        alleleSelect = [alleleSelect ','];
    end
    al = allelelist{i};
    alleleSelect = [alleleSelect al];
    % alleleSelect = 'HLA-DRB1*10:01,HLA-DRB1*12:01';
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

fileID = fopen('IEDBresult.dat','w');
nbytes = fprintf(fileID,'%c',response);
fclose(fileID);

res = readtable('IEDBresult.dat','Delimiter','\t');
for i = 1:nallele
    temp = res(strcmp(res.allele,allelelist{i}),'percentile_rank');
    rank(i,1) = temp.percentile_rank;
end

if nallele == 1
    rank(2,1) = rank(1,1); %Homozygote check
end

% convert rank to affinity
% kon = f(rank);