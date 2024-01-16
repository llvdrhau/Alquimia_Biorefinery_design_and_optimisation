function R =my_kinetics(Opt,R, t)
%% function my_kinetics
% R =my_kinetics(Opt,R, t)
% returns structure R updating:
% - reaction rates (metabolic + active transport)
% - pmf-related ATP rate
% Inputs:   Opt         additional parameter structure
%           t           time
%           R           parameter structure
% Outputs:
%           R           parameter structure

%% Loading variables
Kr = R.pKt.Kin.KrV;    z_r = Opt.z_r;
St = R.St;             pOp = R.pOp;
rm = R.rm;             AlgSt = R.AlgSt;
stoM = full(R.rm.stoM); 
K = R.pKt.Kin.K;

if Opt.electronSource
    eSource=Opt.eSource;
else
    eSource=0;
end

f = AlgSt.f.f_r; 
f = f*St.X;                                                                % (mol X/Lr)
R.rm.f = f;

%% 
R = f_my_regulation(R,Opt,t);                                                  % Determines rate of Anabolism, Decay and Na+-pump
h=R.rm.h;
f=R.rm.f;

%% Loading Monod terms

Monod = Opt.Monod;

% pos_Glu=Opt.pos_Glu;
% Monod(pos_Glu)=(St.Ce_Glu/(1e-3+St.Ce_Glu));


%% Reaction rate determination
r = f.*K.*Kr.*h.*Monod;                                                    % f is the feasibility vector
                                                                           % K is an activation vector (set in excel)
                                                                           % Kr is the maximum rate of the reaction (first term of Monod equation)
                                                                           % h determines the rate of some reactions (Anabolism, Decay and Na_pump)
                                                                           % Monod is the concentration-dependent regulation part of the Monod equation
                           
r = r.*(1 + z_r);                                                          % Application of the reaction optimised vector
r = r*(pOp.Vr/pOp.Vx);                                                     % (mol_i/Lx h).
R.rm.r = r; %update

pos_EMP = find(strcmp(rm.rmNames,'EMP'));
r(pos_EMP) = my_Gly();                                                     %Calculates glycolysis rate
R.rm.r = r;
%% Checking if the system fulfills NADH consrevation
% if abs(sum(r.*stoM(strcmp(St.StNames, 'Ci_NADH'),:)')+eSource*r(pos_EMP)*2) > 1e-3 && t ~= 0 
%     balance_NADH=abs(sum(r.*stoM(strcmp(St.StNames, 'Ci_NADH'),:)')+eSource*r(pos_EMP)*2);    % (mol NAD/Lx h)
%     fprintf('NADH balance is: %f \n',balance_NADH)
%     fprintf('At time: %f h.\n',t)
% end
R.rm.r = r;

%% Active transport
R=f_my_Act_trans(R, Opt);                                                  % Determines the rate of active transport and substrate transport

%% Energetics
indexATP = strcmp(rm.rmNames,'ATPsynthase');
r(indexATP) = f_my_energetics(R);                                          % Determines the pmf-related and maintenance ATP rate

R.rm.r = r;

%% my_Gly
    function fGly = my_Gly()                                               % Determines glycolysis rate
        r = R.rm.r;
        auxGly = find(strcmp(rm.rmNames, 'NADH > H2')); 
        fGly = sum(r(1:auxGly).*-stoM(strcmp(St.StNames, 'Ci_Pyr'),1:auxGly)')/2; % Sum of the pyruvate-consumption reaction rates divided by 2. (mol Glu/Lx h)
    end
end