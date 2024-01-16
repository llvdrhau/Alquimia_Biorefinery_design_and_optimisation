function [der, R, Opt] = my_model(t, x, R)
%% function mymodel
% [sys, R] = my_model(t, x, idR, R)
% returns the derivatives of the states (x) for each time step
% Inputs:   t           time
%           x           states at time t
%           R           structure with all data??
% Outputs:
%           sys         dx/dt  derivatives at time t
%           R           parameters structure R
%           Opt         structure Opt with additional parameters
%% Making again x
act_states = R.St.act_states;
% pos=[1:19 38:46 87:90];
% vec=ones(93,1);
% vec(pos)=0;
states(act_states==1)=x;
states=states';
states(act_states==0)=R.St.StV(act_states==0);
%% R update
if R.flagOpt == 0
    St = R.St;
    pOp = R.pOp;
    iNae= strcmp(St.StNames,'Ce_Na');                                          % Position of external Na+
    iNai= strcmp(St.StNames,'Ci_Na');                                          % Position of internal Na+
                                                                                                                                                                                                                                                          % Sets external Na+ concentration equal to the internal one.
    auxCC = (states <= pOp.Clim);                                                   % States lower than the concentration limit
    states = auxCC.*pOp.Clim + (1-auxCC).*states;                                        % States with lower-than-the-limit concentrations are set to the limit
    StV = states;                                                                   % Update StV with the states values
    
    aux = St.StNames;
    for i=1:(St.numSt)
        St.(char(aux(i))) = states(i);                                              % Update St with the states values
    end
    
    StVchr = [StV(1:St.pos_CitC)                                               % Intracellular transportable states
        1                                                                      % Internal Water concentration
        StV((St.pos_CitC+1):St.pos_Ci)                                         % Intracellular non-transportable states
        StV(St.pos_S:St.pos_Ce)                                                % Biomass and extracellular states
        1                                                                      % External water concentration
        StV((St.pos_Ce+1):end)];                                               % Gas states
    
    %Update in R
    R.St = St;
    R.St.StV = StV;
    R.St.StVchr = StVchr;
end
%% Physicochemical calculations
if R.flagOpt
    R = my_kineticsOpt(R.Opt, R);                                          % For using in optimization mode to calculate the energetics of each indiviual reaction.
end
if R.flagOpt == 0
    R = my_feeding(t, R);                                                  % Updates feeding-related variables
    Inf_ = R.Inf_;                                                         % Influent features
    Xinf = Inf_.InfV(3:end);                                               % Influent concentrations
    
    Q_eff = Inf_.Q;                                                        % Inlet and outlet flows are equal (CSTR)
    
    R = my_algebraics(R,t);                                                % Calculation of deltaG for metabolic and transport reactions and feasiblity factor. Also speciation matrixes.
    
    % Determination of biomass and liquid volumes
    Xt = sum(R.St.StV(strcmp((R.St.Phase),'S')));                          % Biomass concentration (C-mol/L)
    pOp.Xt = Xt;
    Vx = (pOp.Xt*pOp.Vr)/pOp.rho;                                          % Biomass volume (L)
    pOp.Vx = Vx;
    pOp.Vliq = pOp.Vr - Vx;                                                % Liquid volume (L)
    
    if pOp.Vliq <= 0
        error('Error in the input values.-  The operational values introduced are not consistent. Volume of liquid calculated <= 0')
    elseif pOp.Xt <= 0
        error('Error.- Microorganisms concentration is 0')
    end
    
    R.pOp = pOp;                                                           % Update of pOp in R.
    R = my_transport(R);                                                   % Calculation of passive product transport
    %% Optimization section
    Opt = evalin('base', 'Opt');
    %if t >= Opt.T0
    
    my_z(R, t);
    
    %end
    
    %% Mass balances and derivates calculation
    % For each state variable the mass balance is applied according to the variable phase
    
    R = my_kinetics(evalin('base', 'Opt'), R, t);                          % Calculation of reaction, active transport rates and pmf-related energetics
    
    R = my_stoichiometry(R);                                               % Calculation of reaction rates of each compounds
    
    NRt = R.NRt;                                                           % Transport-related compound rates
    NR = R.NR;                                                             % Reaction-related compound rates
    
    St = R.St;
    pOp = R.pOp;
    dX_dt = zeros(St.numSt,1);
    
    % For biomass
    index = strcmp(St.Phase,'S');                                          % Biomass position
    dX_dt(index) = -St.StV(index)/pOp.SRT + NR(index);                     % mol X/Lr·h
    
    %Volume changes are not considered
    pOp.dXt_dt = 0; %TODO Poderase quitar?
    R.pOp = pOp;
    
    % For gaseous compounds. Gas mixing balance in the head space to calculate Qgas accounting for water pressure in equilibrium
    
    NinT = Inf_.Qgas * pOp.Pgas / (pOp.Rg*pOp.T);                          % Moles of gas in the influent (mol_Gas/h)
    ngas_tr = sum(strcmp(St.Phase, 'tG').* NRt) * pOp.Vgas;                % Moles of gas generated (mol_Gas/h)
    
    NoutT = NinT + ngas_tr;                                                % Moles of gas to get out (mol_Gas/h)
    
    Qgas = pOp.Rg * pOp.T * NoutT*(NoutT > 0) / (pOp.Pgas-pOp.Ph2o);       % Calculation of the gas flow (L_Gas/h)
    R.Qgas = Qgas;
    
    % Transportable intracellular compounds
    index = strcmp(St.Phase,'tC');
    dX_dt(index) = -(St.StV(index)/pOp.Xt)*pOp.dXt_dt + NR(index) + NRt(index); % (mol_i/Lx h)
    
    % Non-transportable intracellular compounds
    index = strcmp(St.Phase,'mC');
    dX_dt(index) = -(St.StV(index)/pOp.Xt)*pOp.dXt_dt + NR(index);         % (mol_i/Lx h)
    %     if Opt.electronSource                                                  % Correct mass balance in case of electrofermentation
    %         index_NAD=find(strcmp(St.StNames,'Ci_NAD'));
    %         index_EMP=find(strcmp(R.rm.rmNames,'EMP'));
    %         dX_dt(index_NAD)=dX_dt(index_NAD)-R.rm.r(index_EMP)*2*Opt.eSource;
    %         dX_dt(index_NAD+1)=dX_dt(index_NAD+1)+R.rm.r(index_EMP)*2*Opt.eSource;
    %     end
    
    % Extracellular compounds
    index = strcmp(St.Phase,'tR');
    dX_dt(index) = (1/pOp.Vliq)*(Q_eff*(Xinf(index) - St.StV(index))) + ((St.StV(index)*pOp.Vr)/(pOp.Vliq*pOp.rho))*pOp.dXt_dt + NR(index) + NRt(index); % (mol_i/Lliq h)
    
    % Gas compounds
    index = strcmp(St.Phase,'tG');
    dX_dt(index) = (1/pOp.Vgas)*(Inf_.Qgas*Xinf(index) - Qgas*St.StV(index)) + NRt(index); % (mol_i/Lgas h)
    
    % Inert gas compounds
    index = strcmp(St.Phase,'tGi');
    dX_dt(index) = (1/pOp.Vgas)*(Inf_.Qgas*Xinf(index) - Qgas*St.StV(index)); % (mol_i/Lgas h)
    
    %Derivatives determination
    sys = dX_dt +  auxCC.*(dX_dt < 0).*(-dX_dt);                           % Those compounds under the concentration limit and with negative derivates have a derivative equal to 0
    R.dX_dt = sys;
    der=sys(act_states==1);
    R.t=t;
    assignin('base', 'R', R)
    
else
    der = zeros(length(states),1);                                              % In optimization mode derivates are artificially set to 0.
    R.dX_dt = der;
end


