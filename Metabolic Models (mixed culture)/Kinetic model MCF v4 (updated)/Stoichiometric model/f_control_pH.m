function [Sh, spcM, Ce_Cat, Ce_An] = f_control_pH(R)
%% function f_control_pH
%[Sh, spcM, Ce_Cat, Ce_An] = f_control_pH(R)
% Returns:
%          - extracellular proton concentration
%          - extracellular speciation matrix 
%          - extracellular cation/anion concentration needed to keep extracellular pH constant
% Inputs:   
%           R           parameter structure
% Outputs:
%           Sh          extracelluar proton concentration
%           spcM        extracellular speciation matrix
%           Ce_Cat      extracellular cations concentration
%           Ce_An       extracellular anions concnetration
% Rebeca Gonzalez-Cabaleiro. University of Santiago de Compostela. Spain
% Modified by Alberte Regueira
% Last modification 12/04/2018
%Please contact alberte.regueira@usc.es if you intend to use this code.

%% Loading varaibles
St  = R.St;    StVchr = St.StVchr;
pTh = R.pTh;   Keq = pTh.Keq;
pOp = R.pOp;
StV  = StVchr((St.pos_Schr+1):St.pos_Cechr);      % Extracellular concentrations 
chrM = pTh.chrM(1:length(StV),:);                 % Extracelular charges matrix 
Keq = Keq(1:length(StV),:);                       
w = St.Ce_H2O;

%% Checking the existence of a zero pool in the function between pH 1 and 14
Sh = 10^(-pOp.pH);                                % Extracellular proton concentration
spcM = zeros(size(chrM));
Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);

spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
spcM(:,2) = (StV * Sh^3)                                    ./Denm;
spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
F = Sh + sum(sum(spcM.*chrM));

if F < 0 %If global charge is negative Ce_cat is equal to that value to compensate
    Ce_Cat = abs(F);
    Ce_An = 0;
    
else %If global charge is positive Ce_An is equal to that value to compensate
    Ce_An = F;
    Ce_Cat = 0;
end
end

