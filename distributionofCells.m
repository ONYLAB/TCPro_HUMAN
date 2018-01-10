function NumCellsper_ml = distributionofCells(MinNumPBMCs,MaxNumPBMCs)

bookname = 'PBMCdistribution.xlsx';
opts = detectImportOptions(bookname);
% opts.DataRange='B3';
% opts.VariableNamesRange = 'A2';
opts.RowNamesRange='A3:A8';

CellTypeDist = readtable(bookname,opts);

fun = @(x)(sum(x)-100)^2;

x0 = [-1,2];
lb = CellTypeDist{:,2};
ub = CellTypeDist{:,3};
x0 = rand(6,1).*(ub-lb) + lb; %Random initial Start Point
options = optimoptions('fmincon','Display','off');

x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
% validate x
sum(x);
% 5th row is CD8+ so it will be depleted
totcellafterdepletion = sum(x([1:4 6]));
scale = 100.0/totcellafterdepletion;
x = x*scale;
x(5) = 0;

rng('shuffle');
InitPBMC = rand*(MaxNumPBMCs-MinNumPBMCs) + MinNumPBMCs; %Random initial number of cells/mililiter

NumCellsper_ml = InitPBMC*x/100;