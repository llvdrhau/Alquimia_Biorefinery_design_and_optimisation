function [R,error_lable] = my_algebraics(R,t)
%% function my_algebraics
% R = my_algebraics(idR, R)
% returns structure R updating:
% -the value of pH
% -the speciation matrix R.AlgSt.spcM_IN R.AlgSt.spcM_EX
% -the Gibbs formation matrix
% -feasibility factor
% etc.
% Inputs:   t           time
%           R           parameter structure
% Outputs:
%           R           parameter structure

persistent Sh_inicio  %Kept in memory for function my_algebraics. To be used as an initial guess for the pH calculation
%% Loading variables
St = R.St;    StVchr = St.StVchr;
pTh = R.pTh;   pOp = R.pOp;
stoMfull = R.rm.stoMfull;   % Coefficient matrix for metabolic reactions.
trpMfull = R.rmTr.trpMfull; % Coefficient matrix for cellular transport reactions.
K_active = R.pKt.Kin.K;     % Reactions activated in Excel.

StV_IN = StVchr(1:St.pos_Cichr);                 % Intracelular states
StV_EX = StVchr((St.pos_Schr+1):St.pos_Cechr);   % Extracelular states
chrM_IN = pTh.chrM(1:St.pos_Cichr,:);            % Charges matrix of intracelular states
chrM_EX = pTh.chrM(1:length(StV_EX),:);          % Charges matrix of extracellular states

%% Intracellular pH determination
if size(Sh_inicio) ~= 0
    R.St.Ci_H=Sh_inicio;           % Initial guess from last iteration
end
[Sh_IN, spcM_IN] = f_solve_pH(R);  % Intracellular pH and speciation matrix determination

% if error_Lable == true 
%     errorLable = true ;
% end 

Sh_inicio = Sh_IN;

%% Extracellular pH control
[Sh_EX, spcM_EX, Ce_Cat, Ce_An] = f_control_pH(R);        % Determines the needed acid or base concentrion to keep extracellular pH constant.
St.Ce_Cat = (Ce_Cat ~= 0)*Ce_Cat + (Ce_Cat == 0)*St.Ce_Cat; % If base is added extracelluar cations concentration is updated.
St.StV(strcmp(St.StNames,'Ce_Cat')) = St.Ce_Cat;             % Update Ce_Na in St.StV
St.StVchr(strcmp(St.StchrNames,'Ce_Cat')) = St.Ce_Cat;      % Update Ce_Na in St.StVchr
St.Ce_An = (Ce_An ~= 0)*Ce_An + (Ce_An == 0)*St.Ce_An;    % If acid is added extracellular anions concentration is updated.
St.StV(strcmp(St.StNames,'Ce_An')) = St.Ce_An;            % Update Ce_An in St.StV
St.StVchr(strcmp(St.StchrNames,'Ce_An')) = St.Ce_An;      % Update Ce_An in St.StVchr

%% Determining speciation matrixes

% Without water or protons
spcM = [spcM_IN(1:St.pos_CitC,:);                                                                                          % Intracellular transportable compounds
    spcM_IN((St.pos_CitC+2):St.pos_Cichr,:);                                                                            % Intracellular non-transportable compounds
    [zeros((St.pos_Schr-St.pos_Cichr),1) StVchr((St.pos_Cichr+1):St.pos_Schr) zeros((St.pos_Schr-St.pos_Cichr),3)];     % Biomass
    spcM_EX(1:end-1,:);                                                                                                 % Extracellular compounds
    [zeros(((St.numSt+2)-St.pos_Cechr),1) StVchr((St.pos_Cechr+1):(St.numSt+2)) zeros(((St.numSt+2)-St.pos_Cechr),3)]]; % Gas compounds


% With water and protons
spcMfull = [spcM_IN;
    [0 Sh_IN 0 0 0];                                                                                                    % Intracelllular proton concentration
    [zeros((St.pos_Schr-St.pos_Cichr),1) StVchr((St.pos_Cichr+1):St.pos_Schr) zeros((St.pos_Schr-St.pos_Cichr),3)] ;    % Biomass
    spcM_EX;
    [0 Sh_EX 0 0 0];                                                                                                    % Extracellular proton concentration
    [zeros(((St.numSt+2)-St.pos_Cechr),1) StVchr((St.pos_Cechr+1):(St.numSt+2)) zeros(((St.numSt+2)-St.pos_Cechr),3)] ];% Gas compounds

chrM_all = [pTh.chrM(1:St.pos_Cichr,:);                  % Charges of intracelluar compounds
    [0 1 0 0 0];                                  % Proton charge
    pTh.chrM((St.pos_Cichr+1):St.pos_Schr,:);     % Biomass charge
    pTh.chrM(1:(St.pos_Cechr-St.pos_Schr),:);     % Charges of extracellular compounds
    [0 1 0 0 0];                                  % Proton charge
    pTh.chrM((St.pos_Schr+1):end,:)];             % Gas compounds charges


% Calculation of the ionic strength Ic
Ic_IN = 0.5*(sum(sum([spcM_IN;[0 Sh_IN 0 0 0]].*([chrM_IN;[0 1 0 0 0]].^2))));      % (mol/L)
Ic_EX = 0.5*(sum(sum([spcM_EX;[0 Sh_EX 0 0 0]].*([chrM_EX;[0 1 0 0 0]].^2))));      % (mol/L)

% Calculation of the cell wall osmotic pressure
pi_IN = sum(StV_IN)*pOp.Rg*pOp.T;           % Osmotic pressure in atm (Van't Hoff equation)
pi_OUT = sum(StV_EX)*pOp.Rg*pOp.T;          % Osmotic pressure in atm (Van't Hoff equation)
Pcell = pi_IN - pi_OUT;

%% Gibbs energy calculation for metabolic reactions
spcRfull = accumarray([(1:(St.numSt+4))',pTh.spcRfull],1,[(St.numSt+4),5]);  % Vector containing which species of each compounds is considered to react (with water and protons)
spcR = accumarray([(1:St.numSt)', pTh.spcR],1,[St.numSt,5]);                 % Vector containing which species of each compounds is considered to react
v_r = sum(spcR.*spcM,2);                                                     % Concentration of reacting species of each compound
v_rt = v_r;
% Gibbs formation energy matrix
GfM_all = [pTh.GfM(1:St.pos_Cichr,:);                % Intracellular compounds
    [0 0 0 0 0];                               % Proton
    pTh.GfM((St.pos_Cichr+1):St.pos_Schr,:);   % Biomass
    pTh.GfM(1:(St.pos_Cechr-St.pos_Schr),:);   % Extracellular compounds
    [0 0 0 0 0];                               % Proton
    pTh.GfM((St.pos_Schr+1):end,:)];           % Gas compounds

% DGr0 (standard gibbs energy) with Debye-Hückel correction

A=stoMfull'*sum(spcRfull.*GfM_all,2);   % Dgr0
B=2.303*pOp.Rth*pOp.T*pOp.A;
C = ((sqrt(Ic_IN)/(1+pOp.B*sqrt(Ic_IN)))*(stoMfull(1:St.pos_Cichr+1,:)'*(sum(spcRfull(1:St.pos_Cichr+1,:).*chrM_all(1:St.pos_Cichr+1,:),2)).^2));             % Intracellular charge balance
D = ((sqrt(Ic_EX)/(1+pOp.B*sqrt(Ic_EX)))*(stoMfull((St.pos_Cichr+2):end,:)'*(sum(spcRfull((St.pos_Cichr+2):end,:).*chrM_all((St.pos_Cichr+2):end,:),2)).^2)); % Extracellular charge balance

DGr0 = A-B*(C+D);                                                          % (kJ/mol)

% Correction for real activities

E=pOp.Rth*pOp.T*(stoMfull'*nansum(spcRfull.*log(spcMfull),2));    % Correction for real concentration
F=stoMfull(1:St.pos_Cichr+1,:)'*sum(spcRfull(1:St.pos_Cichr+1,:).*chrM_all(1:St.pos_Cichr+1,:),2)*pOp.F*pOp.E; % Correction for electric field.
DGr = DGr0 + E + F;

DGr(end) = DGr(end) + pOp.pmf;    % Adding a proton extrusion to Na+ inwards transport (Na+ is modelled to be transported antiporter with protons)

spcN = accumarray([(1:St.numSt)', pTh.spcN],1,[St.numSt,5]);               % Neutral species of each compound
v_rtFree = sum(spcN.*spcM,2);                                              % Concentration of the neutral species


%% Gibbs energy calculation for transport reactions

% DGtr0 (standard gibbs energy) with Debye-Hückel correction
G=trpMfull'*sum(spcRfull.*GfM_all,2); % Dgtr0
H=- 2.303*pOp.Rth*pOp.T*pOp.A;
I=(sqrt(Ic_IN)/(1+pOp.B*sqrt(Ic_IN)))*(trpMfull(1:St.pos_Cichr+1,:)'*(sum(spcRfull(1:St.pos_Cichr+1,:).*chrM_all(1:St.pos_Cichr+1,:),2)).^2);             % Intracellular charge balance
J=(sqrt(Ic_EX)/(1+pOp.B*sqrt(Ic_EX)))*(trpMfull((St.pos_Cichr+2):end,:)'*(sum(spcRfull((St.pos_Cichr+2):end,:).*chrM_all((St.pos_Cichr+2):end,:),2)).^2); % Extracellular charge balance

DGtr0=G+H*(I+J);

% Correction for real activities
K=pOp.Rth*pOp.T*(trpMfull'*nansum(spcRfull.*log(spcMfull),2));             % Correction for real concentration
L=trpMfull(1:St.pos_Cichr+1,:)'*sum(spcRfull(1:St.pos_Cichr+1,:).*chrM_all(1:St.pos_Cichr+1,:),2)*pOp.F*pOp.E; % Correction for electric field.

DGtr= DGtr0 + K + L;                                                       % (kJ/mol)

%% Na pump

Antiport = trpMfull(1:St.pos_Cichr+1,:)'*sum(spcRfull(1:St.pos_Cichr+1,:).*chrM_all(1:St.pos_Cichr+1,:),2); % Charge balance for transport reactions
DGtr =  DGtr + Antiport.*pOp.pmf;                                                                           % Transport is modelled as an antirporter with protons.

%% Updating variables
nan_locations = isnan(DGr);
DGr(nan_locations) = 0;
AlgSt.DGr0 = DGr0;
AlgSt.DGr = DGr;
AlgSt.DGtr = DGtr;

%% Feasibility factors
f_r=(tanh(-50*(DGr-pOp.DGmin))+1)/2;                                       % Step function depending on the value of DGr

Opt = evalin('base', 'Opt');
pos_EMP=find(strcmp(R.rm.rmNames,'EMP'));
pos_NADH=find(strcmp(R.rm.rmNames,'NADH > H2'));
K_active=K_active(pos_EMP+1:pos_NADH-1);                                   % Select only catabolic active reactions

A = Opt.f_r(pos_EMP+1:pos_NADH-1);                                         % Feasibility factor of the previous time step
B = f_r(pos_EMP+1:pos_NADH-1);                                             % Feasibility factor of the current time step
C = Opt.z_r(pos_EMP+1:pos_NADH-1)+1;                                       % Reactions selected by the optimizer
if norm((A-B).*(K_active.*C))>1e-3                                         % If the feasibility factor has changed for those reactions selected by the optimizer
    Opt.T0 = t;                                                             % set optimizing time to the current time
end
Opt.f_r=f_r;
assignin('base', 'Opt', Opt)

f_rt=tanh(-10*(DGtr-pOp.DGmin));                                           % Step function depending on the value of DGtr
%% Updating variables

St.Ci_H  = Sh_IN;
St.Ce_H = Sh_EX;

AlgSt.Ic_IN = Ic_IN;
AlgSt.Ic_EX = Ic_EX;
AlgSt.Pcell = Pcell;
AlgSt.spcM_IN = spcM_IN;
AlgSt.spcM_EX = spcM_EX;
AlgSt.spcM = spcM;
AlgSt.v_r = v_r;
AlgSt.v_rt = v_rt;
AlgSt.v_rtFree = v_rtFree;
AlgSt.f.f_r = f_r;
AlgSt.f.f_rt = f_rt;

R.AlgSt = AlgSt;
R.St = St;
R.pOp = pOp;

end








