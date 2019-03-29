function [DC,NK,BC,CD4,CD8,MC] = collectdonorPBMC_analyze(MinNumPBMCs,MaxNumPBMCs,Vp,seed)

Vp = Vp*1e3; %now in ml
NumCellsper_ml(:,1) = distributionofCells(MinNumPBMCs,MaxNumPBMCs,seed);

DC = Vp*NumCellsper_ml(1);
NK = Vp*NumCellsper_ml(2);
BC = Vp*NumCellsper_ml(3);
CD4 = Vp*NumCellsper_ml(4);
CD8 = Vp*NumCellsper_ml(5);
MC = Vp*NumCellsper_ml(6);

% Write out number if desired
% if ~exist('InitNumCell.txt')
%     PBMCcells{1}='DC: ';
%     PBMCcells{2}='NK: ';
%     PBMCcells{3}='BC: ';
%     PBMCcells{4}='CD4: ';
%     PBMCcells{5}='CD8: ';
%     PBMCcells{6}='MC: ';
%     fileID = fopen('InitNumCell.txt','w');
%     fprintf(fileID,'%s\n','Init Cell Counts');
%     for i = 1:6
%         fprintf(fileID,'%s %d\n', PBMCcells{i}, Vp*NumCellsper_ml(i));
%     end
%     fclose(fileID);
% end

% disp('Init Cell Counts');
% disp(['DC: ' num2str(round(DC))]);
% disp(['NK: ' num2str(round(NK))]);
% disp(['BC: ' num2str(round(BC))]);
% disp(['CD4: ' num2str(round(CD4))]);
% disp(['CD8: ' num2str(round(CD8))]);
% disp(['MC: ' num2str(round(MC))]);