function [Response,pval] = Main_human(donor_ID,cohort,SampleConcentration,Fp,pdMS0) 

% Fixed ODESet options
%Make sure that all solution components are nonnegative (13 elements all >=0)
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NonNegative', 1:13);

% Time points for ODE set
tspan1 = 0:0.25:8;

%% AssayType == 'ELISpot'
    
Va = 2e-4; %L
ELISpot = zeros(6,2);

for n = 1:6

    for SimType = 0:1
        
        % Load the parameters
        Parameters(SimType,cohort,Va,donor_ID,SampleConcentration,Fp,pdMS0); %SimType=1 if with Sample, 0 if without
        load Parameters.mat; %#ok<LOAD>
        
        % Initial condition vector
        yic1 = [Ag0; MS0; ID0; MD0; AgE0; pE0; ME0; pME0; pM0; M0; NT0; AT_N0; Prolif0];
        
        % Call ODE
        [T1,Y1]=ode15s(@f, tspan1, yic1, options, pars); 
        
        IL2S = Y1(:,(20+13*N):(23+13*N));
        numIL2secretors = sum(IL2S(end,[1 3:4])); %[MD*RhoNT; DeltaNT*D_N.*Fp.*NT; abs(RhoAT.*E_N.*AT_N)]
        
        ELISpot(n,SimType+1) = numIL2secretors;
        
    end
end
%% AssayType == 'incorporation'
    
%This is the scale when you take aliquots 2ml to 100 mul
Va = 2e-3; %L
VolAliquot = 1e-4; %L
scaleVOL = VolAliquot/Va;

IncubationTime = 18/24; %days

Incorporation = zeros(3,2,4);

for n=1:3

    for SimType = 0:1
        % Load the parameters
        Parameters(SimType,cohort,Va,donor_ID,SampleConcentration,Fp,pdMS0); %SimType=1 if with Sample, 0 if without
        load Parameters.mat; %#ok<LOAD>
        
        % Initial condition vector
        yic1 = [Ag0; MS0; ID0; MD0; AgE0; pE0; ME0; pME0; pM0; M0; NT0; AT_N0; Prolif0];
        
        % Call ODE
        [T1,Y1]=ode15s(@f, tspan1, yic1, options, pars); 
        
        for DayLimit = 5:8
            Incorporation(n,SimType+1,DayLimit-4) = getincubation(DayLimit,IncubationTime,Y1,T1,scaleVOL,Va,options,pars);
        end
    end

end

for i = 1:4
    sample = squeeze(Incorporation(:,2,i));
    medium = squeeze(Incorporation(:,1,i));
    Response(i) = mean(sample)/mean(medium); %#ok<AGROW>
    [~,pval(i)] = ttest2(medium,sample,'Tail','left','Vartype','unequal'); %#ok<AGROW>
end

sample = ELISpot(:,2);
medium = ELISpot(:,1);
Response(5) = mean(sample)/mean(medium);
[~,pval(5)] = ttest2(medium,sample,'Tail','left','Vartype','unequal');


function IncorporationResponse = getincubation(DayLimit,IncubationTime,y_record,t_record,scaleVOL,Vp,options,pars)

load('Parameters.mat','N'); %Comes from ELISPOT run, we already looked at this there

numberoftimesamples = 10;

[~,t_index] = min((t_record-DayLimit).^2);

% Get Initial conditions
Ag1 = y_record(t_index,1);
MS1  = y_record(t_index,2); %No new Maturation signal
ID1 = y_record(t_index,3);
MD1 = y_record(t_index,4);
AgE1 = y_record(t_index,5);
pE1 = y_record(t_index,(5+1):(5+N));
ME1 = y_record(t_index,(6+N):(11+N));
pME1 = y_record(t_index,(12+N):(11+7*N));
pM1 = y_record(t_index,(12+7*N):(11+13*N));
M1 = y_record(t_index,(12+13*N):(17+13*N));
NT1 = y_record(t_index,18+13*N);
AT_N1 = y_record(t_index,(19+13*N):(18+14*N));
Prolif1 = y_record(t_index,(19+14*N):(19+17*N));

% Initial condition vector
Prolif1 = 0.0*Prolif1; %This is a place holder for proliferation
yic2 = [Ag1, MS1, ID1, MD1, AgE1/scaleVOL, pE1/scaleVOL,  ME1/scaleVOL, pME1/scaleVOL,  pM1/scaleVOL, M1/scaleVOL,NT1, AT_N1, Prolif1]';
yic2 = yic2*scaleVOL;

%pars, Vp changes
pars(3)=Vp*scaleVOL;

% dosing interval
t_start = DayLimit;
t_end = DayLimit+IncubationTime;
tspan2 = linspace(t_start, t_end, numberoftimesamples);

% Call ODE
[T2,Y2] = ode15s(@f, tspan2, yic2, options,pars); %#ok<ASGLU>

% Prolif, activated helper T cells derived from memory T cells, number of cells
Prolif_vector = Y2(:,(20+13*N):(21+13*N)); %MD*RhoNT; sAt.*RhoAT.*E_N.*AT_N;

IncorporationResponse = sum(Prolif_vector(end,:)); %#Proliferated cells