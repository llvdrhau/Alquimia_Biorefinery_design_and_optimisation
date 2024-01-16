function [Sh, spcM] = f_solve_pH(R)
%% function f_solve_pH
%[Sh, spcM] = f_solve_pH(R)
% Returns:
%          - intracellular proton concentration
%          - intracellular speciation matrix 
% Inputs:   
%           R           parameter structure
% Outputs:
%           Sh          proton concentration
%           spcM        intracellular speciation matrix
% Rebeca Gonzalez-Cabaleiro. University of Santiago de Compostela. Spain
% Modified by Alberte Regueira
% Last modification 12/04/2018
%Please contact alberte.regueira@usc.es if you intend to use this code.

%% Loading variables
St  = R.St;      StVchr = St.StVchr;
pTh = R.pTh;     Keq = pTh.Keq;
Sh_inicio = St.Ci_H;
StV  = StVchr(1:St.pos_Cichr);       % Intracelular concentrations
chrM = pTh.chrM(1:St.pos_Cichr,:);   % Intracelular charges matrix 
w = St.Ci_H2O;                       % Water concentration

%% Checking the existence of a zero pool in the function between pH 1 and 14

a=1e-14;
b=1;
Sh_v = [a; b]; F = zeros(2,1);
for i = 1:2
    Sh=Sh_v(i);
    spcM = zeros(size(chrM));
    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
    
    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm; %Not hydrated form concentrations     
    spcM(:,2) = (StV * Sh^3)                                    ./Denm; %Fully protonated form concentrations
    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm; %1st deprotonated form concentrations
    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm; %2nd deprotonated form concentrations
    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm; %3rd deprotonated form concentrations
    
    F(i) = Sh + sum(sum(spcM.*chrM));                                   % Evaluation of the charge balance for the current Sh value, F(Sh)
end

if prod(F) > 0
    error_lable = true;  %gooood help me plz
    %error('ERROR.- The sum of charges returns a wrong value')
end

fa = F(1);
fb = F(2);

%% pH determination (Newton-Raphson method)
Sh = Sh_inicio;
i=1;    Tol = 5.e-15;    maxIter = 50;
spcM = zeros(size(chrM));      % Inicialization of matrix of species
dspcM = zeros(size(chrM));
while i <= maxIter
    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
    spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
    
    
    F = Sh + sum(sum(spcM.*chrM));  % Evaluation of the charge balance for the current Sh value, F(Sh)
    
    % Calculation of all derivated functions
    
    dDenm = Denm.^2;
    aux = 3*Sh^2*(Keq(:,1)/w + 1) + 2*Sh*Keq(:,2) + Keq(:,2).*Keq(:,3);
    
    dspcM(:,1) =  (3*Sh^2*Keq(:,1).*StV)./(w*Denm) - ((Keq(:,1).*StV*Sh^3).*aux) ./(w*dDenm);
    dspcM(:,2) =  (3*Sh^2*StV)./Denm - (StV*Sh^3.*aux) ./dDenm;
    dspcM(:,3) = (2*Sh*Keq(:,2).*StV)./Denm - ((Keq(:,2).*StV*Sh^2).*aux) ./dDenm;
    dspcM(:,4) = (Keq(:,2).*Keq(:,3).*StV)./Denm - ((Keq(:,2).*Keq(:,3).*StV*Sh).*aux)./dDenm;
    dspcM(:,5) = -(Keq(:,2).*Keq(:,3).*Keq(:,4).*StV.*aux) ./dDenm;
    
    dF = 1 + sum(sum(dspcM.*chrM));     % Evaluation of the charge balance for the current Sh value, dF(Sh)
    err = F/dF;                         % Error
    Sh = Sh - err;                      % Newton-Raphson algorithm
    
    if (abs(err) < 1e-14) && (abs(F) < Tol)
        if (Sh > 1e-14) && (Sh < 1)    % Checking if a valid pH was obtained
            return
        else
            break
        end
    end
    i = i+1;
end

%% pH determination (Pegasus method)

i = 1; maxIter = 50;
n1 = 0; n2 = 0;
while (i < maxIter)
    Sh = (fb*a-fa*b)/(fb-fa);
    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
    
    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
    spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
    
    fc = Sh + sum(sum(spcM.*chrM));
    if fa*fc > 0
        n1 = n1+1;
        if n1 == 2
            fb = (fc/(fc+fa))*fb;
            n1 = 0;
        end
        a = Sh; fa = fc;
    elseif fb*fc > 0 % To avoid problems when fc == 0
        n2 = n2+1;
        if n2 == 2
            fa = (fc/(fc+fb))*fa;
            n2 = 0;
        end
        b = Sh; fb = fc;
    end
    
    err1 = abs(fc);
    err2 = abs(Sh-(fb*a-fa*b)/(fb-fa));
    if (err1 < Tol) && (err2 < 1e-14)
        return
    end
    i = i+1;
end
error('ERROR.- The interaction loop for pH calculation has not got convergence')
end