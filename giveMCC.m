function [MCC,cm] = giveMCC(Lamberth,Drug)
% Give Matthews Correlation Coefficient

[c,cm,ind,per] = confusion(Lamberth,Drug);

FN = cm(2,1);
FP = cm(1,2);
TN = cm(1,1);
TP = cm(2,2);

MCC = TP*TN - FP*FN;
a = TP + FP;
b = TP + FN;
c = TN + FP;
d = TN + FN;

MCC = MCC/sqrt(a*b*c*d);

if isnan(MCC) 
    MCC=-1; 
end

