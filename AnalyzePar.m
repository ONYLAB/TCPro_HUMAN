seqtable = readtable('runpeptides.dat');
changes = [0.5 2];
changingparameters();
parnames = readtable('changingparameters.txt');

for k = 1:2
    
    ParChange = changes(k);
    
    for j = 1:length(parchangelist)
        
        ParChangeIndex = parchangelist{j};
        
        for i = 1:height(seqtable)
            
            seqid = seqtable{i,'SEQNAME'}{1};
            epitope = seqtable{i,'SEQUENCE'}{1};
            disp(seqid);
            mkdir(seqid)
            cd(seqid)
            [responsesummary,significantresponsesummary,kon,responsevector,significancevector] = runalldonor(epitope,ParChangeIndex,ParChange);
            kons{i,j,k}=kon;
            allmeanresponses{i,j,k} = responsevector;
            allsignificance{i,j,k} = significancevector;
            CUMresponses(i,j,k,:)=responsesummary;
            CUMsignificance(i,j,k,:)=significantresponsesummary;
            cd ..
            save
            
        end
    end
end

for k = 1:2
    
    ParChange = 1;
    
    for j = 1:1
        
        ParChangeIndex = parchangelist{j};
        
        for i = 1:height(seqtable)
            
            seqid = seqtable{i,'SEQNAME'}{1};
            epitope = seqtable{i,'SEQUENCE'}{1};
            disp(seqid);
            mkdir(seqid)
            cd(seqid)
            [responsesummary,significantresponsesummary,kon,responsevector,significancevector] = runalldonor(epitope,ParChangeIndex,ParChange);
            allmeanresponses0{i,j,k} = responsevector;
            allsignificance0{i,j,k} = significancevector;
            CUMresponses0(i,j,k,:)=responsesummary;
            CUMsignificance0(i,j,k,:)=significantresponsesummary;
            cd ..
            
        end
    end
end

close all
for i = 1:height(seqtable)
    
    for j = 1:length(parchangelist)
        
        for k = 1:2
            
            responses(i,j,k) = allmeanresponses{i,j,k}{1}(4);%Day 8 of the radiactive proliferation experiment
            responses0(i,j,k) = allmeanresponses0{i,1,k}{1}(4);%Day 8 of the radiactive proliferation experiment
            difference(i,j,k) = 100*(responses(i,j,k)-responses0(i,j,k))/responses0(i,j,k);
            
        end
        
        sortchange(i,j) = sum(abs(difference(i,j,:)));
        
    end
end

for i = 1:height(seqtable)
    
    figure
    [~,Ind] = sort(sortchange(i,:),'descend');
    
    for j = 1:length(parchangelist)
        
        for k = 1:2
            
            if k ==1
                if j ==1
                    scatter(j,difference(i,Ind(j),k),50,'b','filled');
                else
                    scatter(j,difference(i,Ind(j),k),50,'b','filled','HandleVisibility','off');
                end
                hold on
            else
                if j ==1
                    scatter(j,difference(i,Ind(j),k),50,'r','filled');
                else
                    scatter(j,difference(i,Ind(j),k),50,'r','filled','HandleVisibility','off');                    
                end
                hold on
            end
            
        end
        
        if j ==1
            legend('50% of Original Parameter Value','200% of Original Parameter Value','Location','NorthEast');            
        end
        
    end
    xticks(1:length(parchangelist));
    xticklabels(parnames.Params(Ind));
    xtickangle(45)
    LT = 3; %Line thickness
    AxFS = 24; %Ax Fontsize
    AxLW = 2; %Ax LineWidth
    xlabeltext = 'Parameters';
    ylabeltext = '%Stimulation Index Change';
    set(gcf,'color','w');
    set(gca,'fontsize', AxFS);
    set(gca,'LineWidth',AxLW);
    xlabel(xlabeltext);
    ylabel(ylabeltext);
    box on
    xlim([0 length(parchangelist)+1])
    plot(0:length(parchangelist)+1,zeros(1,length(parchangelist)+2),'k--','HandleVisibility','off')
    title(['Donor 1 on Antigen Peptide: ' seqtable{i,'SEQNAME'}{1}])
end