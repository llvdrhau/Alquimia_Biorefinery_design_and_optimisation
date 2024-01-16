function [F_acidogenic_fermentation,End_concentrations_M] = f_determine_my_F(R, Opt, pH, consumed_AA, consumed_glucose)


% This function determines the acidogenic fermentation stoichiometric 
% parameters in the stationary state
% *********************************************************************** %
% INPUTS:
% D
% pH
% VECTOR of all AA in mol/l !!
% glucose in mol/l!!
% *********************************************************************** %
% OUTPUT:
% Acidogenic fermentation compounds in mol/l !!  
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


%% Modifying the feeding

D = [1/24];      % Dilution rate (h-1)

%% Modifying model parameters

R.pOp.Gan=[5e3];             % Gain of the active transport control (f_my_Act_trans.m)
R.pOp.Ana=5;           % Gain of the glucose anabolism control (f_my_regulation.m)
R.pOp.control=[1e-3];        % Tolerance of the optimiser to do a new reaction selection process (my_z.mat)

%% Modifying simulation parameters

R.St.StV(find(strcmp('Ce_Arg',R.St.StNames)):find(strcmp('X',R.St.StNames))) = 1e-3;  % AA intracelular concentration (M)
R.St.StV(strcmp('X',R.St.StNames)) = 1e-2;    % Biomass initial concentration (M)
R.St.StV(strcmp('Ce_Glu',R.St.StNames)) = 1e-19;  % External glucose concentration (M). Set to 1e-19 when simulation a proten monofermentation, otherise set to 1e-4.


%% ODE options

x0 = R.St.StV;
vec = R.St.act_states;
x0=x0(vec==1);
nneg = 1:length(x0);
tspan = 0:5:150;
options = odeset('OutputFcn', @(t, x, flag) out(t, x, flag), 'RelTol', R.TolR,'AbsTol', R.TolA,'NonNegative', nneg, 'MaxStep', 0.01); %options for the ode solver
clc
%sound(sin(1:3000));
fprintf('\n> Bioenergetic model is running >>>>>')

% A version of the original R and Opt are saved to use at the begginig of each loop
R_ori=R;      
Opt_ori=Opt;

%% Simulation

tic
delete All_StatesVar.mat
R=R_ori;
Opt=Opt_ori;

% Feeding update
pos1 = strmatch('Ce_Arg',R.St.StNames)+2;
posend = strmatch('Ce_Hist',R.St.StNames)+2;
R.Inf_.InfV(pos1:posend)=consumed_AA; 
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
R.pOp.pH=pH;
R.pOp.pmf=0.2*96.5+298.15*8.314e-3*log(10^(-pH)/1e-7);

[T,Y] = ode15s(@(t, x) my_model(t, x, R), tspan, x0, options);  % sates are in mol/L 
            
%sound(sin(1:3000));


%% Stoichiometric acidogenic fermentation factors calculation

load('All_StatesVar')

Qgas=All_StatesVar(end,384);
Qliq=R.Eff.Q;

nSt=R.St.numSt;
y=All_StatesVar(:,2:nSt+1);
Ac=y(:,strcmp(R.St.StNames,'Ce_Ac'));
Pro=y(:,strcmp(R.St.StNames,'Ce_Pro'));
But=y(:,strcmp(R.St.StNames,'Ce_Bu'))+y(:,strcmp(R.St.StNames,'Ce_iBu'));
Val=y(:,strcmp(R.St.StNames,'Ce_Val'))+y(:,strcmp(R.St.StNames,'Ce_iVal'));
iCap=y(:,strcmp(R.St.StNames,'Ce_iCap'));
Et=y(:,strcmp(R.St.StNames,'Ce_EtOH'));
H2=y(:,strcmp(R.St.StNames,'Ce_For'))+y(:,strcmp(R.St.StNames,'Ce_H2'))+y(:,strcmp(R.St.StNames,'G_H2'))*Qgas/Qliq;



aux=Ac(end)+Pro(end)+But(end)+Val(end)+iCap(end)+Et(end)+H2(end);

leg_names={'Acetate','Propionate','Butyrate','Valerate','H2','Ethanol'};

F_matrix={};
F_matrix(:,1)=leg_names;
F_matrix(2:end+1,:) = F_matrix(1:end,:);
F_matrix(1,:) = {''};

F_matrix{1,2}='End OUT Conc mol/l (M)';
F_matrix{2,2}=Ac(end);
F_matrix{3,2}=Pro(end);
F_matrix{4,2}=But(end);
F_matrix{5,2}=Val(end);
F_matrix{6,2}=H2(end);
F_matrix{7,2}=Et(end);

F_matrix{1,3}='End OUT molar fraction';
F_matrix{2,3}=Ac(end)/(aux-iCap(end));
F_matrix{3,3}=Pro(end)/(aux-iCap(end));
F_matrix{4,3}=But(end)/(aux-iCap(end));
F_matrix{5,3}=Val(end)/(aux-iCap(end));
F_matrix{6,3}=H2(end)/(aux-iCap(end));
F_matrix{7,3}=Et(end)/(aux-iCap(end));

COD=[96,112,160,208,16,64]; %arriba de todo
COD=COD';
aux=[F_matrix{2:end,3}]';
aux=aux.*COD;

F_matrix{1,4}='End OUT COD fraction';
F_matrix{2,4}=aux(1)/sum(aux);
F_matrix{3,4}=aux(2)/sum(aux);
F_matrix{4,4}=aux(3)/sum(aux);
F_matrix{5,4}=aux(4)/sum(aux);
F_matrix{6,4}=aux(5)/sum(aux);
F_matrix{7,4}=aux(6)/sum(aux);

F_acidogenic_fermentation=[F_matrix{2:7,4}]';  % COD fractiom

biomass = y(:,strcmp(R.St.StNames,'X')); 
biomassEnd = biomass(end); 

End_concentrations_M = [F_matrix{2:7,2},biomassEnd]';  % mol/L
end

