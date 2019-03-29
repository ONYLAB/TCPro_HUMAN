function NumCellsper_ml = distributionofCells(MinNumPBMCs,MaxNumPBMCs,seed)

fun = @(x)(sum(x)-100)^2; %x is the vector of cell percentages

%Donor cell distribution variability fixed
rng(seed);

% Dendritic Cells ; NK Cells; B Cells; CD4+; CD8+; Monocytes
lb = [0.63;6.84;6.36;26.15;12.38;6.04];
ub = [1.46;35.17;16.51;48.26;33.11;11.56];
x0 = rand(6,1).*(ub-lb) + lb; %Random initial Start Point fixed for cohort member

options = optimoptions('fmincon','Display','off');

%Variation in the cell type distribution
x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
% validate x, make sure it is 100% total
while sum(x)<99.9 || sum(x)>100.1    
    x = fmincon(fun,x,[],[],[],[],lb,ub,[],options);
end
    
% 5th row is CD8+ so it will be depleted
totcellafterdepletion = sum(x([1:4 6]));
scale = 100.0/totcellafterdepletion;
x = x*scale;
x(5) = 0;

%Variation in the assay
rng('shuffle');
rng('shuffle');
InitPBMC = rand*(MaxNumPBMCs-MinNumPBMCs) + MinNumPBMCs; %Random initial number of cells/mililiter

NumCellsper_ml = InitPBMC*x/100;