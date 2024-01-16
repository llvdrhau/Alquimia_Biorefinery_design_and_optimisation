function R =f_my_regulation(R,Opt,t)
%% function my_regulation
% This function 
% h= my_regulation(R, idR)
%  updates vector h(r???)  and f(Favourability multiplied by biomass conc.)
% Calculates r(ana), r(decay) and r(Na_pump)


%% Loading parameters
rm = R.rm;  %Parameters related to metabolic reactions (coefficints matrix, etc)
St = R.St; %States
AlgSt = R.AlgSt;
f = AlgSt.f.f_r; %Vector of favourable reactions.
f = f*St.X; 

%% GLYCOLYSIS REGULATION.- Equal to Glucose transport rates WHY IS IT ZERO? 

iEMP = strcmp(rm.rmNames,'EMP'); %3rd reaction  WHY IS IT ZERO? I guess it is because the EMP ratio is already fixed.
f(iEMP) = 0; %Set to 0 element number 3 in f

%% Na+ pump
iNa = strcmp(rm.rmNames, 'Na_Pump'); %Position of reaction Na_Pump
Na_pump = R.rm.h(iNa); %Inizialization value. It is only updated at the beginning of the simulation. Afterwards it is allways the same.
h = ones(rm.num_r,1); %Column vector of ones
h(iNa) = Na_pump; %Vector columns of 1s and h value for Na_pump at the end
f(iNa) = St.X; %It is forced to be positive. WHY? I guess it is because it has a already fixed ratio.

%Proportional control action of the Na_Pump. See Section G.
% if St.Ci_H == 1e-7
%     qNa = 0; 
% else
%     aux = St.Ci_H > 1e-7; %aux=1 if pH<7. aux=0 if pH>7.
%     %qNa = aux*(1-exp(K_.K_Na*((1e-7-St.Ci_H)))) + (1-aux)*(exp(K_.K_Na*((St.Ci_H-1e-7)))-1); %Value of qNa+ depending on pH medium. Proportional control.
%     qNa = (1-aux)*(exp(K_.K_Na*(1e-7-1e-14/St.Ci_H))-1) + aux*(1-exp(K_.K_Na*(1e-7-St.Ci_H)));
% end
pH = -log10(St.Ci_H);
Gan = 1e-2*tanh(10*(7-pH));
qNa = Gan*(pH-7)^2;
h(iNa) =  h(iNa) + qNa; %mol Na/mol Cx·h

% Corections for h(iNa)
DGrATP = R.AlgSt.DGr(strcmp(rm.rmNames,'ATPsynthase')); %deltaG of ATP formation (UNITS??)

if DGrATP < 40 && sign(h(iNa)) ~= AlgSt.f.f_r(iNa) %If deltaG(ATP)<40 (there is a lot of available ADP) and the sign of h(Na_Pump) and favorability of Na_pump is different then
    h(iNa) = h(iNa)*(DGrATP^10/(DGrATP^10 + 15^10)); %Correction in case DGrATP is very low (<20 kJ/mol). It slows down the Na_pump to save energy (to confirm)
end

if DGrATP > 80
    h(iNa) = h(iNa)*(1-DGrATP^20/(DGrATP^20 + 80^20));
end

if h(iNa) < 0 %If Na+ is being pumped out
    h(iNa) = h(iNa)*(St.Ci_Na^2/(St.Ci_Na^2 + 1e-12)); %Regulation. If no internal Na, the reaction is slowed down
else          %If Na+ is being pumped in
    h(iNa) = h(iNa)*(St.Ce_Na^2/(St.Ce_Na^2 + 1e-12)); %Regulation. If no external Na, the reaction is slowed down
end

if h(iNa)>10          %For extreme cases (big pH differences) it increases Na_pump to fix the situation quicker. 
    h(iNa)=h(iNa)*10; %mol Na/mol Cx·h
elseif h(iNa)<-10
    h(iNa)=h(iNa)*(-10); %mol Na/mol Cx·h
end

%% ANABOLISM GLUCOSE

iAnab = strcmp(rm.rmNames, 'Anab');
f(iAnab) = St.X;
ana = (DGrATP>=50)*((DGrATP-50)/R.pOp.Ana); %Equation [S45]
ana = ana*Opt.Monod_term(1); %Equation [S46] 
ana = min([ana 1]);  % Maximum anabolism is 1 (minimum value between 1 and the calculated value)
h(iAnab) = ana;  %Equation [S46] %mol Cx/mol Cx·h

%% ANABOLISM PROTEINS

iAnabP = strcmp(rm.rmNames, 'AnabProt');
f(iAnabP) = St.X;
ana = (DGrATP>=50)*((DGrATP-50)/5); %Equation [S45]
pos1=find(strcmp(St.StNames,'Ce_Arg'));
posend=find(strcmp(St.StNames,'Ce_Hist'));
nC=[6 3 4 6 5 3 4 3 2 5 5 6 6 5 5 4 6];
C_prot=mean(St.StV(pos1:posend)'.*nC);
ana = ana*(C_prot/(6e-3+C_prot)); %Equation [S46] 
ana = min([ana 1]);  % Maximum anabolism is 1 (minimum value between 1 and the calculated value)
h(iAnabP) = ana;  %Equation [S46] %mol Cx/mol Cx·h



%% DECAY Glu
iDec = strcmp(rm.rmNames,'DecayGlu');
f(iDec) = St.X;
decay = (DGrATP<50)*((50-DGrATP)/R.pOp.Ana); %If DgrATP is higher than 50 kJ/mol there is no decay. %Equation [S45]
h(iDec) = min([decay 1]); %mol ??/mol Cx·h
if DGrATP < 0
    fprintf('DGrATP very low!! \n')
    fprintf('At time: %f h.\n',t)
end
%% DECAY Prot
iDec = strcmp(rm.rmNames,'DecayProt');
f(iDec) = St.X;
decay = (DGrATP<50)*((50-DGrATP)/5); %If DgrATP is higher than 50 kJ/mol there is no decay. %Equation [S45]
h(iDec) = min([decay 1]); %mol ??/mol Cx·h


%% Update

R.rm.h = h;  %Returns vector h. All values are 1 except for positions 1 (Ana) 2 (decay) and 40 (Na_pump)
R.rm.f = f;  %Updated f for Na, EMP, anabolism and decay

