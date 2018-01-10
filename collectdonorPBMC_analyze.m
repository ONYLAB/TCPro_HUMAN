function [DC,NK,BC,CD4,CD8,MC] = collectdonorPBMC_analyze(MinNumPBMCs,MaxNumPBMCs,Vp)

Vp = Vp*1e3; %now in ml
NumCellsper_ml(:,1) = distributionofCells(MinNumPBMCs,MaxNumPBMCs);

DC = Vp*NumCellsper_ml(1);
NK = Vp*NumCellsper_ml(2);
BC = Vp*NumCellsper_ml(3);
CD4 = Vp*NumCellsper_ml(4);
CD8 = Vp*NumCellsper_ml(5);
MC = Vp*NumCellsper_ml(6);

% disp('Init Cell Counts');
% disp(['DC: ' num2str(round(DC))]);
% disp(['NK: ' num2str(round(NK))]);
% disp(['BC: ' num2str(round(BC))]);
% disp(['CD4: ' num2str(round(CD4))]);
% disp(['CD8: ' num2str(round(CD8))]);
% disp(['MC: ' num2str(round(MC))]);