function dydt=f(t,y,pars) %#ok<INUSL>

% Parameter vector
NA=pars(1);
Dose=pars(2); %#ok<NASGU>
Vp=pars(3);
%% T-epitope characteristics of therapeutic proteins
N=pars(4);
kon=reshape(pars(5:(4+6*N)),N,6);
koff= reshape(pars((5+6*N):(4+12*N)),N,6);
%% Dendritic cells
BetaID=pars(5+12*N);
DeltaID=pars(6+12*N);
KMS=pars(7+12*N);
BetaMD=pars(8+12*N);
%% Antigen presentation
AlphaAgE=pars(9+12*N);
BetaAgE=pars(10+12*N);
Beta_p=pars(11+12*N);
BetaM=pars(12+12*N);
Beta_pM=pars(13+12*N);
kext=pars(14+12*N);
kin=pars(15+12*N);
KpM_N=pars(16+12*N);
VD=pars(17+12*N);
VE=pars(18+12*N);
%% T helper cells
RhoNT=pars(19+12*N);
BetaNT=pars(20+12*N);
DeltaNT=pars(21+12*N);
RhoAT=pars(22+12*N);
BetaAT=pars(23+12*N);
Fp=pars(24+12*N);
%% Initial conditions as parameters in the differential equations:
ID0=pars(25+12*N); %#ok<NASGU>
ME0=pars((26+12*N):(31+12*N));
EpitopeRatio=pars(32+12*N);

%% State variables for the differential equations
% Ag, antigenic protein in the plasma compartment
Ag=y(1);
% MS, maturation signal for Immature dendritic cells
MS=y(2);
% ID, Immature dendritic cells
ID=y(3);
% MD, mature dendritic cells
MD=y(4);
%AgE, Ag in the endosomes
AgE=y(5);
% pE, free epitope peptides from Ag digestion (N Epitopes)
pE=y((5+1):(5+N));
% ME(1:6), free MHC II molecule in endosome (6 MHCIIs allowed)
ME=y((6+N):(11+N));
% pME(N,6),	epitope peptide-MHC II complex in endosomes
pME=reshape(y((12+N):(11+7*N)),N,6); % pME is in a matrix form, N epitopes against 6 possible MHC alleles
% pM(N,6), epitope peptide-MHC II complexes on dendritic cell membrane
pM=reshape(y((12+7*N):(11+13*N)),N,6);
% M, free MHC II molecule on dendritic cell membrane (6 MHCIIs allowed)
M=y((12+13*N):(17+13*N));
% NT, na�ve helper T cells
NT=y(18+13*N);
% AT_N activated helper T cells (N Epitopes)
AT_N=y((19+13*N));
% Prolif, Placeholder for proliferating cells MD+AT
Prolif=y((20+13*N):(23+13*N)); %#ok<NASGU>

%% Calculate functions for helper T cells activation, proliferation, or differentiation

% Adding up all the pM molecules from one epitope against 6 different MHC alleles
pM_NUMBER=EpitopeRatio*sum(pM*ones(6,1)*1E-12*NA); %pM is an Nx6 matrix times 6x1 gives Nx1 summing it gives 1x1 scalar number

% The activation function D for helper T cells
D_N=(MD/(MD+sum(Fp.*NT)+sum(AT_N)))*(((pM_NUMBER)./(pM_NUMBER+KpM_N)));

% The proliferation/differentiation function E for helper T cells
E_N=(MD/(MD+sum(Fp.*NT)+sum(AT_N)))*((pM_NUMBER-KpM_N)./(pM_NUMBER+KpM_N));

%% Differential equations
% Ag, y(1), total amount of antigenic protein in the well, pmole
dydt(1,1)=-(ID+MD)*AlphaAgE*VD*(Ag/Vp);

% MS, y(2), maturation signal, particularly, LPS, for immature dendritic cells, ng
dydt(2,1)=-(ID+MD)*AlphaAgE*VD*(MS/Vp);

% ID, y(3), immature dendritic cells, cells
dydt(3,1)=-BetaID*ID-DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS);

% MD,	y(4), mature dendritic cells, cells
dydt(4,1)=DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS)-BetaMD*MD;

%AgE, y(5), Ag in the endosome, pmole
dydt(5,1)=AlphaAgE*VD*(Ag/Vp)-BetaAgE*AgE;

% pE, y((5+1):(5+N)), T-epitope peptides from Ag digestion, pmole
dydt((5+1):(5+N),1)=BetaAgE*AgE-Beta_p*pE-kon.*(pE*(ME'/VE))*ones(6,1)+koff.*pME*ones(6,1);

% ME, y((6+N):(11+N)),free MHC-II molecule in endosome, pmole
dydt((6+N):(11+N),1)= BetaM*(ME0-ME)-(ones(1,N)*(kon.*(pE*(ME'/VE))))'+(ones(1,N)*(koff.*pME))'+kin*M;

% pME, y((12+N):(11+7*N)), T-epitope-MHC-II complex in endosome, pmole
dydt((12+N):(11+7*N),1)= reshape(kon.*(pE*(ME'/VE))-koff.*pME -Beta_pM*pME -kext*pME,6*N,1);

% pM,	y((12+7*N):(11+13*N)), T-epitope-MHC-II complex on dendritic cell membrane, pmole
dydt((12+7*N):(11+13*N),1)=reshape(kext*pME -koff.*pM,6*N,1);

% M, y((29+13*N):(34+13*N)), free MHC-II molecule on dendritic cell membrane, pmole
dydt((12+13*N):(17+13*N),1)=-kin*M +(ones(1,N)*(koff.*pM))';

% NT, y(18+13*N), na�ve helper T cells pool, cells
dydt(18+13*N,1)=MD*RhoNT-BetaNT*NT-sum(DeltaNT*D_N.*Fp.*NT);

% AT_N, y((19+13*N):(18+14*N)), activated helper T cells derived from NT cells
dydt((19+13*N),1)=DeltaNT*D_N.*Fp.*NT+RhoAT.*E_N.*AT_N-BetaAT*AT_N; 

% Prolif, placeholder for proliferating cells MD+AT, y((19+14*N):(18+15*N)) 
sAt = E_N>=0; %Making sure E_N is positive so we take into account only the proliferating cells
dydt((20+13*N):(23+13*N),1)=[MD*RhoNT; sAt.*RhoAT.*E_N.*AT_N; DeltaNT*D_N.*Fp.*NT; abs(RhoAT.*E_N.*AT_N)];