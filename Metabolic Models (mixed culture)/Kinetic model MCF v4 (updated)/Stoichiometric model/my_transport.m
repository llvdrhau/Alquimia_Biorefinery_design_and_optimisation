function R = my_transport(R)
%% function my_transport
% R = my_transport(R)
% Returns:
%          - passive transport rates
% Inputs:   
%           R           parameter structure
% Outputs:
%           R           parameter structure
% Rebeca Gonzalez-Cabaleiro. University of Santiago de Compostela. Spain
% Modified by Alberte Regueira
% Last modification 12/04/2018
%Please contact alberte.regueira@usc.es if you intend to use this code.

%% Loading variables
pOp = R.pOp;         St = R.St;       pTh = R.pTh;
trpM = R.rmTr.trpM; 
spcM = R.AlgSt.spcM;
spcR = R.pTh.spcR;
K_dif = R.pKt.KTr.K_difV;          %Diffusivity constants for liquid transport reactions (L/molCx h) and kla (h-1) for gas transport.
StTrans = R.AlgSt.v_rtFree;        %Concentration of neutral species of the compounds

aux = strcmp(St.StNames,'Ci_CO2') + strcmp(St.StNames,'Ce_CO2'); 
StTrans = (1-aux).*StTrans + aux.*(spcM(:,1)+spcM(:,2));                   % For CO2 the hydrated form is added to the neutral form.
% StTrans(find(aux)) = spcM(find(aux),1);

w = strcmp(St.Phase, 'tG');                                                % Gas compounds from reactions                                     
StGT = zeros(St.numSt,1);

%Correction for CO2 liquid-gas transport
spcR(strcmp(St.StNames,'Ce_CO2')) = 1;                                     % Change the reactive species of extracellular CO2 to hydrated form.

for i=1:St.numSt                                                           % StGas contains the liquid concentrations of compounds being transferred to the gas phase
    if w(i) ~= 0 
        namegas = St.StNames{i};
        namegas = namegas(3:end);
        name =horzcat('Ce_',namegas); 
        j = strcmp(St.StNames,name); 
        u = j.*spcM(:,spcR(j)); 
        StGT = StGT + u; 
    end
end

ntG = sum(w);                                                              % Number of liquid-gas transport reactions
StGT = StGT + w.*(pTh.Kh.KhV.*StTrans*pOp.Rg*pOp.T);                       % Add to StGT gas concentraions (Henry's Law) (Calculate the liquid concentrations in eq with the gas)
aux = [zeros((length(K_dif)-ntG),1);ones(ntG,1)]; 
rt = (1-aux).*(((-trpM')*StTrans).*K_dif*St.X) + aux.*(K_dif.*((-trpM')*StGT));   % transport rates calculation (mol/Lr h for liquid-liquid transport and mol/Lliq·h for liquid-gas transport)

%% Variables update
R.rmTr.rt_dif = rt;