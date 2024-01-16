function [z] = f_my_optimization(R, Opt, t0, RUTA, TRANSP)
%% function my_optimization
% z = my_optimization(idR, R, Opt, t0, RUTA, TRANSP)
% returns the reaction optimisation vector
% Inputs:   R           parameter structure
%           Opt         additional parameter structure
%           RUTA        Matrix for change basis between the vectorial space of all reactions and the space of eligible in optimisation reactions
%           TRNASP      Matrix linking which transport reactions are activated by each of the reactions eligible in optimisation
% Outputs:
%           z           Reaction optimisation vector

%% Initializing reactions
St = R.St;
R.flagOpt = 1; 
RUTAtot = RUTA;
TRANSPtot = TRANSP;
index_nz=Opt.index_nz;
x0 = St.StV;
act_states = R.St.act_states;
% pos=[1:19 38:46 87:90];
% vec=ones(93,1);
% vec(pos)=0;
x0=x0(act_states==1);
dATP = zeros(1, size(RUTAtot,1)); 
dNADH = zeros(1, size(RUTAtot,1)); 
%dCO2 = zeros(1, size(RUTAtot,1)); 
R.Opt = Opt;
RR = struct(R);
nz = size(RUTAtot,1); 

%% Generating the optimising vectors
for i = 1: nz                                                              % A structure RR is created for each eligible reaction assuming that all the substrate 
    RR(i) = R;                                                             % (glucose or an individual AA) is consumed in each of the reactions
    RUTA = RUTAtot(i,:)'; 
    TRANSP = TRANSPtot(i,:)'; 
    RR(i).Opt.z_r=RUTA-1; 
    RR(i).Opt.z_rt=TRANSP-1; 
end
f_r=RUTAtot*R.AlgSt.f.f_r;

pos1 = find(strcmp(St.StNames,'Ce_Glu'));
posend = find(strcmp(St.StNames,'Ce_Hist'));
feeding = St.StV([pos1 pos1+2:posend]);
no_feeding = find(feeding<1e-8);
% no_feeding(length(no_feeding)+1) = 3;  % Deactivate Alanine consumption
desact = ones(length(dATP),1);
for i=1:length(no_feeding)
    if no_feeding(i)==1
        f_r(1:index_nz(1)) = 0;
        desact(1:index_nz(1)) = 0;
    else
        f_r(index_nz(no_feeding(i)-1)+1:index_nz(no_feeding(i))) = 0;
        desact(index_nz(no_feeding(i)-1)+1:index_nz(no_feeding(i))) = 0;
    end
end
optimise = find(desact==1);

for k = 1:length(optimise)
    j=optimise(k);                                                            % For each of the eligible reactions run my_model in optimisation mode to calculate the energetics
    [~, RS] = my_model(t0, x0, RR(j));                                     %  of that singular reaction 
    
    if Opt.electronSource                                                  % In case of electro-fermentation
        eSource = Opt.eSource;
    else
        eSource=0;
    end
    a = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_NADH'), 1:end).*RS.rm.r(1:end)');              % NADH rate (mol NADH/Lx h)
    b = -(sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_Pyr'), 1:end).*RS.rm.r(1:end)'));            % Pyruvate rate (only active in glucose degradation) (mol Pyr/Lx h)
    
    dATP(j) = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_ATP'), 1:end).*RS.rm.r(1:end)')+b;         % ATP rate (mol ATP/Lx h)
    %dCO2(j) = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_CO2'), 1:end).*RS.rm.r(1:end)');         % CO2 rate (mol CO2/Lx h)
    dNADH(j) = (1+eSource)*b + a;
        
end

dATP((dATP<-100))=-100;              
dNADH(isinf(dNADH)) = -100;
dATP(isnan(dATP)) = 0;
dNADH(isnan(dNADH)) = 0;
%dCO2(isnan(dCO2)) = 0;
f_r(index_nz)=1;                                                    % Force all null reactions to be eligible

%% Optimiser
try
    [z,EXITFLAG,~,~] = f_quadraticProgram(Opt.z,dATP,dNADH,nz,index_nz,f_r);
catch
    [z,EXITFLAG,~] = f_linProgram(Opt.z,dATP,dNADH,index_nz,f_r);       % In case the quadratic optimiser has any trouble, there is a linear optimiser.
end

if EXITFLAG ~= 1
    fprintf('OLLO!! No solution found for optimisation \n') 
    keyboard
end

