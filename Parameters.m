function Parameters(SimType,epitopes,HLA_DR)

%% Celltype distribution
MinNumPBMCs = 4e6;
MaxNumPBMCs = 6e6;
[ID0,NK,BC,CD4,CD8,MC] = collectdonorPBMC_analyze(MinNumPBMCs,MaxNumPBMCs);
% ID0: the initial number of immature dendritic cell

%% Therapeutic protein PK parameters
% Vp: well volume
Vp = 2e-3;%2.75; % L This is the Volume of the sample
SampleConcentration = 0.3e-6; %Molar This is the actual concentration of the samples with respect to the cell-excluded volume
Dose = SimType*SampleConcentration*Vp*1e12;  % pmole

Endotoxin = 0.0042*1e3; %ng/L
% NA: Avogadro constant
NA = 6.0221367e23;

%% T-epitope characteristics of therapeutic proteins
% N: the number of T-epitope (protein-specific)
% ME0:  initial amount of MHC-II molecule in a single mature dendritic cell pmole
% kon: on rate for for T-epitope-MHC-II binding
% koff:	off rate for for T-epitope-MHC-II binding

% Number of HLA molecules
numHLADR = 1e5;

[kon,koff,N,ME0] = doNETMHCIIpan(epitopes,HLA_DR,SimType,NA,numHLADR);

%% Dendritic cells
% BetaMS: decay rate for the maturation signals (LPS)
BetaMS=0.0;  % day-1

% BetaID:	death rate of immature dendritic cells.
BetaID=0.0924; % day-1

% DeltaID:	maximum activation rate of immature dendritic cells
DeltaID=1.66; % day-1

% Km,MS: LPS concentration at which immature DC activation rate is 50% maximum.
KMS= 9.852E3; % ng/L

% BetaMD:	death rate of mature dendritic cells
BetaMD=0.2310; % day-1

%% Antigen presentation

% konN: on rate for competing peptide-MHC-II binding
konN=ones(6,1)*8.64*1E-3; %  pM-1day-1 NOT USED

% koffN: off rate for competing peptide-MHC-II binding
koffN=8.64*1E-3*ones(6,1)*4000; %  day-1 NOT USED

% AlphaAgE:  Ag internalization rate constant by mature dendritic cells
AlphaAgE=14.4; %day-1

%	BetaAgE	: degradation rate for AgE in acidic vesicles
BetaAgE = 17.28; % day-1
if length(epitopes{1})<25
    BetaAgE = BetaAgE*100;
end

%	Beta_p:	degradation rate for T-epitope peptide (same for all peptides)
Beta_p= 14.4; % day-1

% BetaM: degradation rate for MHC-II
BetaM= 1.663; % day-1.

% Beta_pM: degradation rate for pME
Beta_pM= 0.1663; % day-1

% kext: exocytosis rate for peptide-MHC-II complex from endosome to cell menbrane
kext = 28.8;  % day-1

% kin: internalization rate for peptide-MHC-II complex on DC membrane
kin = 14.4;  % day-1

% KpM_N: number of peptide-MHCII to achieve 50% activation rate of na�ve helper T cells
KpM_N = 400;  % number of complex

% VD: volume of a dendritic cell
VD=0.8463e-012; % L

% VE: volume of endosomes in a dendritic cell
VE=0.05*VD; % L


%% T helper cells
% RhoNT: maximum proliferation rate for activated helper T cells
RhoNT=0.1432; % day-1

%BetaNT: death rate of naive helper T cells
BetaNT=0.3048;%%PREVIOUSLY:0.018; % day-1

%DeltaNT: maximum activation rate of naive helper T cells
DeltaNT=1.5; % day-1

% RhoAT: maximum proliferation rate for activated helper T cells
RhoAT=1.5; % day-1

%BetaAT: death rate of activated helper T cells
BetaAT=0.18; % day-1

%% Initial conditions for state variables in the differential equations:
% Ag0: intitial therapeutic protein in the plasma compartment
Ag0=Dose; % pmole

% MS0: Maturation signal (MS), particularly, endotoxin, LPS
MS0=Endotoxin*Vp;%5*70; % ng

% MD0: the initial number of mature dendritic cells
MD0=0; % cells

% cp0: amount of endogenous competing protein in the plasma
cp0=0.0;%110E6*Vp; % pmole

% cpE0: initial amount of endogenous competing protein in endosome
cpE0=0.0;%1e6/NA*1e12;%110E6*2.75;   % pmole

% cptE0: initial amount of endogenous competing peptide in endosome
cptE0=0;   % pmole

% cptME0: initial amount of endogenous peptide-MHC complex in endosome
cptME0=zeros(6,1); % pmole

% cptM0: initial amount of endogenous peptide-MHC complex on dendritic cell membrane
cptM0=zeros(6,1);  % pmole

% AgE0: initial amount of Ag in endosome
AgE0=0; % pmole

% pE0:	initial amount of T-epitope peptides from Ag digestion in endosome
pE0=ones(N,1)*0;  % pmole

% pME0:	initial amount of T-epitope-MHC-II complex in endosome
pME0=ones(6*N,1)*0;  % pmole

% pM0:	initial amount of T-epitope-MHC-II complex on dendritic cell membrane
pM0=ones(6*N,1)*0;  % pmole

% M0, free MHC-II molecule on dendritic cell membrane
M0=ones(6,1)*0; % pmole

% NT0: the initial number of naive T helper cell = total T helper cells * frequency
% of T-epitope-specific T helper cells (protein-specific)
% Human: Total naive T helper cells per kg body weight is: 876E6 cells/L*5L(blood
% volume)= 4.3800e+009 cells
Fp = ones(N,1)*0.1/1E6;
NT0=CD4;%4.3800e+009*0.33/1E6; % cells

% AT_N0: initial number of activated T helper cell derived from naive T cells
AT_N0=ones(N,1)*0.0; % cells

% AT_M0: initial number of activated T helper cell derived from memory T cells
AT_M0=ones(N,1)*0.0; % cells

% MT0: initial number of memory T helper cell
MT0=ones(N,1)*0.0; % cells

% FT0: initial number of functional T helper cell
FT0=ones(N,1)*0.0; % cells

%% Parameter vector
pars(1)=NA;
% Therapeutic protein PK parameters
pars(2)=Dose;
pars(3)=Vp;
% T-epitope characteristics of therapeutic proteins
pars(4)=N;
pars(5:(4+6*N))=reshape(kon,6*N,1);
pars((5+6*N):(4+12*N))= reshape(koff, 6*N,1);
% Dendritic cells
pars(5+12*N)=BetaMS;
pars(6+12*N)=BetaID;
pars(7+12*N)=DeltaID;
pars(8+12*N)=KMS;
pars(9+12*N)=BetaMD;
% Antigen presentation
pars((10+12*N):(15+12*N))=konN;
pars((16+12*N):(21+12*N))=koffN;
pars(22+12*N)=AlphaAgE;
pars(23+12*N)=BetaAgE;
pars(24+12*N)=Beta_p;
pars(25+12*N)=BetaM;
pars(26+12*N)=Beta_pM;
pars(27+12*N)=kext;
pars(28+12*N)=kin;
pars(29+12*N)=KpM_N;
pars(30+12*N)=VD;
pars(31+12*N)=VE;
% T helper cells
pars(32+12*N)=RhoNT;
pars(33+12*N)=BetaNT;
pars(34+12*N)=DeltaNT;
pars(35+12*N)=RhoAT;
pars(36+12*N)=BetaAT;
pars(37+12*N:36+13*N)=Fp;
% Initial conditions as parameters in the differential equations:
pars(37+13*N)=ID0;
pars(38+13*N)=cpE0;
pars((39+13*N):(44+13*N))=ME0;
pars((45+13*N):(44+14*N))=NT0;
pars((45+14*N):(44+15*N))=MT0;
pars=pars';

save Parameters.mat

function [kon,koff,N,ME0] = doNETMHCIIpan(epitopes,HLA_DR,SimType,NA,numHLADR)
%% Collect -on, -off rates, number of epitopes and amount of initial MHC
% molecules
EpitopePresent = 1;
N = length(epitopes);
if N<1
    EpitopePresent = 0;
    epitopes{1} = 'AAAAAAAAAAA'; %Dummy PolyA epitope
    N=1;
end

if strcmp(HLA_DR{1,1},HLA_DR{1,2}) %homozygot
    disp('Homozygote');
    ME0=[numHLADR; 0.0;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
    % kon: on rate for for T-epitope-MHC-II binding
    kon=EpitopePresent*SimType*repmat([1 0 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
    HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
else
    ME0=[numHLADR/2; numHLADR/2;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
    % kon: on rate for for T-epitope-MHC-II binding
    kon=EpitopePresent*SimType*repmat([1 1 0 0 0 0],N,1)*8.64*1E-3; %  pM-1day-1
    HLAtext = [HLA_DR{1,1} ',' HLA_DR{1,2}];
end

Affinity_DPQ=[4000 4000 4000 4000]; %place holder for other alleles (NOT USED)
Affinity_DR = givekon(epitopes{1},HLAtext);
%	koff:	off rate for for T-epitope-MHC-II binding
koff=8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3; %  day-1

for i = 2:N
    Affinity_DR = givekon(epitopes{i},HLAtext);
    temp = 8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3;
    koff = [koff;temp];
end


function Affinity_DR = givekon(epitopesequence,HLAtext)
%% Run NetMHCIIPan if not done before, extract association rates
% data(1).Sequence = 'epitopesequence'
% data(1).Header = 'Seq'
% fastawrite('example.fsa',data);
if exist('out.dat')==0 %#ok<EXIST>
    dlmwrite('example.fsa',epitopesequence,'')
    command = ['netMHCIIpan -f example.fsa -inptype 1 -xls -xlsfile out.dat -a ' HLAtext];
    [s,m] = unix(command);
end
table = readtable('out.dat');
Affinity_DR(1) = table{1,'nM'};
Affinity_DR(2) = table{1,'nM_1'};
