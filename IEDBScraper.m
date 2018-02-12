function IEDBScraper(peptidesequence,nallele,allelelist)

global all15mers DRBalleles thresholds;

alleleSelect = '';
for i = 1:nallele
    globalDRBindex(i) = find(strcmp(DRBalleles,allelelist{i}));
    if i>1
        alleleSelect = [alleleSelect ','];
    end
    allelelist{i} = strrep(['HLA-' allelelist{i}],'_','*');
    allelelist{i} = [allelelist{i}(1:end-2) ':' allelelist{i}(end-1:end)];
    alleleSelect = [alleleSelect allelelist{i}];
    % alleleSelect = 'HLA-DRB1*10:01,HLA-DRB1*12:01';
end

if exist('IEDBresult.dat')==0 %#ok<EXIST>
    url = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/';
    method = 'method';
    methodSelect = 'recommended';
    sequence_text = 'sequence_text';
    sequence_textSelect = peptidesequence;
    allele = 'allele';
    
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
end

res = readtable('IEDBresult.dat','Delimiter','\t');
for i = 1:nallele
    temp = res(strcmp(res.allele,allelelist{i}),'percentile_rank');
    rank(i,1) = mean(temp.percentile_rank);
    [~,ind] = min((thresholds-rank(i,1)).^2);
    if rank(i,1)<=thresholds(ind)
        th(i,1) = thresholds(ind);
        kd(i,1) = all15mers(globalDRBindex(i),ind);
    elseif rank(i,1)>thresholds(ind) && ind<100
        th(i,1) = thresholds(ind+1);
        kd(i,1) = all15mers(globalDRBindex(i),ind+1);
    else
        th(i,1) = 100;
        kd(i,1) = 5e4;
    end
    
    if nallele == 1
        rank(2,1) = rank(1,1); %Homozygote check
        th(2,1) = th(1,1);
        kd(2,1) = kd(1,1);
    end
end

% convert rank to affinity and write to file
dlmwrite('IEDB_KD.dat',kd');
% kon = f(rank);
% 'IEDB_KD.dat' dlmwrite