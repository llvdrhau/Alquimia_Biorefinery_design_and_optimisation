
load('R Generation/R.mat')
load('R Generation/Opt.mat')

%% Modifying the feeding

% Reading AA profile from excel file
% [AAprofile] = xlsread('Aminoacidos.xlsx','AA'); % Reading AA profile from excel file
% Prot=AAprofile(:,3)';
% conc_glu = [0.018771];
% conc_prot = [10]; % g/L
D = [1/24];      % Dilution rate (h-1)

%% Modifying model parameters

Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
R.pOp.Gan=[5e3];             % Gain of the active transport control (f_my_Act_trans.m)
R.pOp.Ana=5;           % Gain of the glucose anabolism control (f_my_regulation.m)
R.pOp.control=[1e-3];        % Tolerance of the optimiser to do a new reaction selection process (my_z.mat)
RET = [0.8];           % Fraction of ATP yielded in the branched AA reactions (valine, isoleucine, leucine).

%% Modifying simulation parameters

pH =[7];
R.St.StV(find(strcmp('Ce_Arg',R.St.StNames)):find(strcmp('X',R.St.StNames))) = 1e-3;  % AA intracelular concentration (M)
R.St.StV(strcmp('X',R.St.StNames)) = 1e-2;    % Biomass initial concentration (M)
R.St.StV(strcmp('Ce_Glu',R.St.StNames)) = 1e-19;  % External glucose concentration (M). Set to 1e-19 when simulation a proten monofermentation, otherise set to 1e-4.

%% Modyfing stoichiometric matrix
% RET modification in branched AA

%The ATP, ADP, Pi, H2O and H+ stoichiometry is modified in the stoichiometry matrix depending on the RET value. 
% ATTENTION. It is manually adjustes to reactions 1 and 2 of valine, isoleucine and leucine. If additional reactions or some reactions are
% erased the allocation of the variables should be modified


pos_H2O = strcmp(R.St.StNamesFull,'Ci_H2O');
pos_Pi = find(strcmp(R.St.StNamesFull,'Ci_Pi'));
pos_Pi2 = find(strcmp(R.St.StNames,'Ci_Pi'));

pos_Val = strmatch('Val',R.rm.rmNames);
pos_Ile = strmatch('IsoL',R.rm.rmNames);
pos_Leu = strmatch('Leu',R.rm.rmNames);

% Valine
R.rm.stoMfull(pos_H2O,pos_Val(1))=RET-3;  % Water
R.rm.stoMfull(pos_Pi,pos_Val(1))=-RET;    % Pi
R.rm.stoMfull(pos_Pi+1,pos_Val(1))=RET;   % ATP
R.rm.stoMfull(pos_Pi+2,pos_Val(1))=-RET;  % ADP
R.rm.stoMfull(pos_Pi+3,pos_Val(1))=2-RET; % H+

R.rm.stoMfull(pos_H2O,pos_Val(2))=RET-2;  % Water
R.rm.stoMfull(pos_Pi,pos_Val(2))=-RET;    % Pi
R.rm.stoMfull(pos_Pi+1,pos_Val(2))=RET;   % ATP
R.rm.stoMfull(pos_Pi+2,pos_Val(2))=-RET;  % ADP
R.rm.stoMfull(pos_Pi+3,pos_Val(2))=2-RET; % H+

R.rm.stoM(pos_Pi2,pos_Val(1):pos_Val(2))=-RET;  % Pi
R.rm.stoM(pos_Pi2+1,pos_Val(1):pos_Val(2))=RET; % ATP
R.rm.stoM(pos_Pi2+2,pos_Val(1):pos_Val(2))=-RET;% ADP

% Isoleucine
R.rm.stoMfull(pos_H2O,pos_Ile(1))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Ile(1))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Ile(1))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Ile(1))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Ile(1))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Ile(2))=RET-2;
R.rm.stoMfull(pos_Pi,pos_Ile(2))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Ile(2))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Ile(2))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Ile(2))=2-RET;

R.rm.stoM(pos_Pi2,pos_Ile(1):pos_Ile(2))=-RET;
R.rm.stoM(pos_Pi2+1,pos_Ile(1):pos_Ile(2))=RET;
R.rm.stoM(pos_Pi2+2,pos_Ile(1):pos_Ile(2))=-RET;

% Leucine
R.rm.stoMfull(pos_H2O,pos_Leu(1))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Leu(1))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Leu(1))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Leu(1))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Leu(1))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Leu(2))=RET-2;
R.rm.stoMfull(pos_Pi,pos_Leu(2))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Leu(2))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Leu(2))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Leu(2))=2-RET;

R.rm.stoM(pos_Pi2,pos_Leu(1):pos_Leu(2))=-RET;
R.rm.stoM(pos_Pi2+1,pos_Leu(1):pos_Leu(2))=RET;
R.rm.stoM(pos_Pi2+2,pos_Leu(1):pos_Leu(2))=-RET;


%% Initialising ode
x0 = R.St.StV;
vec = R.St.act_states;
x0=x0(vec==1);
nneg = 1:length(x0);
tspan = 0:1/10:100;
options = odeset('OutputFcn', @(t, x, flag) out(t, x, flag), 'RelTol', R.TolR,'AbsTol', R.TolA,'NonNegative', nneg, 'MaxStep', 0.01); %options for the ode solver
clc
fprintf('\n> MODEL RUNNING >>>>>')

% A version of the original R and Opt are saved to use at the begginig of each loop
R_ori=R;      
Opt_ori=Opt;

%% Simulation
% for j=1:1                       % Select manually which columns of AAprofile you want ot use in the simulation
    for i=1:length(pH)
        tic
        delete All_StatesVar.mat
        R=R_ori;
        Opt=Opt_ori;
        % Feeding update       
        pos1 = strmatch('Ce_Arg',R.St.StNames)+2;
        posend = strmatch('Ce_Hist',R.St.StNames)+2;
%         R.Inf_.InfV(pos1:posend)=AAprofile(:,j)/5*conc_prot; % 1= caseína, 2=xelatina
        R.Inf_.InfV(pos1:posend)=consumed_AA; % 1= caseína, 2=xelatina
        pos_Glu = strmatch('Ce_Pyr',R.St.StNames)+2-1;
        R.Inf_.InfV(pos_Glu)=consumed_glucose;
        R.Inf_.InfV(2)=D*2;
        R.pOp.SRT=1/D;
        R.mFd.FdV(2) = D*2;
        R.Eff.EffV(2) = D*2;
        R.Eff.Q = D*2;
        R.Inf_.Q = D*2;
        R.mFd.Q = D*2;
        % pH update
        R.pOp.pH=pH(i);
        R.pOp.pmf=0.2*96.5+298.15*8.314e-3*log(10^(-pH(i))/1e-7);
        
        try
            [T,Y] = ode15s(@(t, x) my_model(t, x, R), tspan, x0, options);
            name=horzcat('All_StatesVar_',num2str(i),'.mat');
            save(name, 'All_StatesVar')
            M = excel(R,Opt);
            name=horzcat('M_',num2str(i),'.mat');
            save(name, 'M')
        catch
            
        end
    end
% end
