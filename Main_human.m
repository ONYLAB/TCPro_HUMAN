% This code were applied for simulating in vivo immune response in human against Adalimumab.
% The dosing regimen was based on clinical dosing regimen of Adalimumab.
function response = Main_human(DayLimit,SimType,epitopes,HLA_DR,donor_ID) %#ok<STOUT>

% clc
close all
% clear all
rng(donor_ID);

% Main simulation

%This is the scale when you take aliquots 2ml to 100 mul
% VolSample = 2000; %muL
% VolAliquot = 100; %muL
% scaleVOL = 1;%VolAliquot/VolSample;

numberoftimespamples = 100;

% Load the parameters
Parameters(SimType,epitopes,HLA_DR); %SimType=1 if with Sample, 0 if without
load Parameters.mat; %#ok<LOAD>

% Run the ODEs
options = odeset('RelTol',1e-10, 'AbsTol',1e-10);
% DayLimit = 8;
IncubationTime = 18/24; %days
Frequency = 1; %Day-1

% First dose
% Initial condition vector
yic1 = [Ag0;MS0; ID0; MD0; cpE0; cptE0;cptME0;cptM0;AgE0;pE0; ME0;pME0; pM0; M0;NT0; AT_N0; AT_M0;MT0; FT0];

% First dosing interval, protein specific
tspan1 =linspace(0,1/Frequency,numberoftimespamples);

% Call ODE
[T1,Y1]=ode15s(@f, tspan1, yic1, options,pars); %#ok<NODEF>

% Record the results into new matrix
t_record=T1;
y_record=Y1;

% Subsequent doses
for num_dose=1:DayLimit  % number of doses, protein specific
    
    % Get Initial conditions
    Ag1=y_record(numel(t_record),1);
    MS1=y_record(numel(t_record),2); %No new Maturation signal
    ID1=y_record(numel(t_record),3);
    MD1=y_record(numel(t_record),4);
    cpE1=y_record(numel(t_record),5);
    cptE1=y_record(numel(t_record),6);
    cptME1=y_record(numel(t_record),(7:12));
    cptM1=y_record(numel(t_record),(13:18));
    AgE1=y_record(numel(t_record),19);
    pE1=y_record(numel(t_record),(19+1):(19+N));
    ME1=y_record(numel(t_record),(20+N):(25+N));
    pME1=y_record(numel(t_record),(26+N):(25+7*N));
    pM1=y_record(numel(t_record),(26+7*N):(25+13*N));
    M1=y_record(numel(t_record),(26+13*N):(31+13*N));
    NT1=y_record(numel(t_record),(32+13*N):(32+13*N));
    AT_N1=y_record(numel(t_record),(33+13*N):(32+14*N));
    AT_M1=y_record(numel(t_record),(33+14*N):(32+15*N));
    MT1=y_record(numel(t_record),(33+15*N):(32+16*N));
    FT1=y_record(numel(t_record),(33+16*N):(32+17*N));
    
    % Initial condition vector
    AT_M1 = 0.0*AT_M1; %This is a place holder for proliferation
    %     MT1 = 0.0*MT1; not 0 now, it's the rest of the naive bystanding t-cells
%     FT1 = 0.0*FT1;
    yic2 = [Ag0, MS0, ID1, MD1, cpE1, cptE1,cptME1,cptM1,AgE1,pE1,  ME1,pME1,  pM1, M1,NT1, AT_N1,AT_M1,MT1,FT1]';
%     yic2 = yic2;
    
    if num_dose ~= DayLimit        
        Duration = 1/Frequency; %Frequency in Days
    else
        Duration = IncubationTime; %IncubationTime in Days        
    end

    % dosing interval
    t_start=num_dose/Frequency;
    t_end=t_start+Duration;
    tspan2 = linspace(t_start, t_end, numberoftimespamples);
    
%     %pars, Vp changes
%     pars(3)=Vp;
%     pars(37+13*N)=ID0;
%     pars((45+13*N):(44+14*N))=NT0;
%     pars((45+14*N):(44+15*N))=MT0;
    
    % Call ODE
    [T2,Y2]=ode15s(@f, tspan2, yic2, options,pars);
    
    % Record the result into new matrix
    t_record=[t_record;T2];
    y_record=[y_record; Y2];
    
end


% Output all state variables
% Ag, Ag in the plasma compartment, pmole
Ag_vector=y_record(:,1); %#ok<NASGU>
% MS, Maturation signal for activation of immature dendritic cells
% (Here is LPS, ng)
MS_vector=y_record(:,2); %#ok<NASGU>
% ID,	Immature dendritic cells, number of cells
ID_vector=y_record(:,3);
% MD,	Mature dendritic cells, number of cells
MD_vector=y_record(:,4);
% cpE, endogenous competing protein in endosome, pmole
cpE_vector=y_record(:,5); %#ok<NASGU>
% cptE, endogenous competing peptide in endosome, pmole
cptE_vector=y_record(:,6); %#ok<NASGU>
% cptME, endogenous peptide-MHC complex in endosome, pmole
cptME_vector=y_record(:,(7:12)); %#ok<NASGU>
% cptM, endogenous peptide-MHC complex on dendritic cell membrane, pmole
cptM_vector=y_record(:,(13:18)); %#ok<NASGU>
%AgE, Ag in the endosomes, pmole
AgE_vector=y_record(:,19); %#ok<NASGU>
% pE,	free epitope peptides from Ag digestion , pmole
pE_vector=y_record(:,(19+1):(19+N)); %#ok<NASGU>
% ME,	free MHC II molecule in endosome , pmole
ME_vector=y_record(:,(20+N):(25+N)); %#ok<NASGU>
% pME,	T-epitope-MHC-II complex in endosome, pmole
pME_vector=y_record(:,(26+N):(25+7*N)); %#ok<NASGU>
% pM,	T-epitope-MHC-II on dendritic cell membrane, pmole
pM_vector=y_record(:,(26+7*N):(25+13*N));
% M, free MHC II molecule on dendritic cell menbrane, pmole
M_vector=y_record(:,(26+13*N):(31+13*N)); %#ok<NASGU>
% NT,	Na�ve helper T cells, number of cells
NT_vector=y_record(:,(32+13*N):(32+13*N));
% AT_N,	activated helper T cells derived from naive T cells, number of cells
AT_N_vector=y_record(:,(33+13*N):(32+14*N));
% AT_M,	activated helper T cells derived from memory T cells, number of cells
AT_M_vector=y_record(:,(33+14*N):(32+15*N));
% MT, memory helper T cells, number of cells
MT_vector=y_record(:,(33+15*N):(32+16*N)); %#ok<NASGU>
% FT, functional helper T cells, number of cells
FT_vector=y_record(:,(33+16*N):(32+17*N)); %#ok<NASGU>

% Post-processing calculation
% Antigen presentation processes

pM_NUMBER_M1=pM_vector(:,(1:N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 1
pM_NUMBER_M2=pM_vector(:,(N+1):(2*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 2
pM_NUMBER_M3=pM_vector(:,(2*N+1):(3*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 3
pM_NUMBER_M4=pM_vector(:,(3*N+1):(4*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 4
pM_NUMBER_M5=pM_vector(:,(4*N+1):(5*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 5
pM_NUMBER_M6=pM_vector(:,(5*N+1):(6*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 6

Total_pM(:,1:N)=pM_NUMBER_M1+pM_NUMBER_M2+pM_NUMBER_M3+pM_NUMBER_M4+pM_NUMBER_M5+pM_NUMBER_M6;   %#ok<NASGU> % Total number of T-epitope-MHC II complex

cptM_NUMBER_M1=cptM_vector(:,(1:N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 1
cptM_NUMBER_M2=cptM_vector(:,(N+1):(2*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 2
cptM_NUMBER_M3=cptM_vector(:,(2*N+1):(3*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 3
cptM_NUMBER_M4=cptM_vector(:,(3*N+1):(4*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 4
cptM_NUMBER_M5=cptM_vector(:,(4*N+1):(5*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 5
cptM_NUMBER_M6=cptM_vector(:,(5*N+1):(6*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 6

Total_cptM(:,1:N)=cptM_NUMBER_M1+cptM_NUMBER_M2+cptM_NUMBER_M3+cptM_NUMBER_M4+cptM_NUMBER_M5+cptM_NUMBER_M6;   %#ok<NASGU> % Total number of T-epitope-MHC II complex


% Save the results
savefile='results.mat';

save(savefile, 'koff', 't_record','Ag_vector','MS_vector','ID_vector','MD_vector',...
    'cpE_vector', 'cptE_vector', 'cptME_vector', 'cptM_vector', 'AgE_vector','pE_vector', ...
    'ME_vector','pME_vector', 'pM_vector', 'M_vector', 'NT_vector', 'AT_N_vector', 'AT_M_vector', 'MT_vector',...
    'FT_vector','Total_pM');

% % % % % % % % % % % % % % % % %
% figure
% LT = 3; %Line thickness
% AxFS = 24; %Ax Fontsize
% AxLW = 2; %Ax LineWidth
% xlabeltext = 'Time (days)';
% ylabeltext = '#cells';
% legenddata = {'Immature DC','Mature DC','Naive T','Activated T'};
% 
% plot(t_record,ID_vector(:,1),'LineWidth',LT);
% hold on
% plot(t_record,MD_vector(:,1),'LineWidth',LT);
% plot(t_record,NT_vector(:,1),'LineWidth',LT);
% plot(t_record,AT_N_vector(:,1),'LineWidth',LT);
% % plot(t_record,MT_vector(:,1),'LineWidth',LT);
% 
% set(gcf,'color','w');
% set(gca,'fontsize', AxFS);
% set(gca,'LineWidth',AxLW);
% legend(legenddata,'Location','eastoutside');
% xlabel(xlabeltext);
% ylabel(ylabeltext);
% set(gca,'yscale','log')
% axis square

response = AT_M_vector(end,1);
numactivatedT = AT_N_vector(end,1) + MD_vector(end,1);