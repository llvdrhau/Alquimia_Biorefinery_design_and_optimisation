function R =my_kineticsOpt(Opt, R)
%% function my_kineticsOpt
% R =my_kineticsOpt(Opt,R)
% returns structure R updating:
% - reaction rates (metabolic + active transport)
% - pmf-related ATP rate
% Inputs:   Opt         additional parameter structure
%           R           parameter structure
% Outputs:
%           R           parameter structure

%% Loading variables
Kr = R.pKt.Kin.KrV;    z_r = Opt.z_r;
St = R.St;             pOp = R.pOp;
rm = R.rm;             AlgSt = R.AlgSt;
K = R.pKt.Kin.K;

f = AlgSt.f.f_r;
f = f*St.X;
R.rm.f = f;

hOpt = my_regulation();                                                    % Sets Anabolism, Decay and Na-pump rates to 0.
Monod = Opt.Monod;

%% Reaction rate determination
r = f.*K.*Kr.*hOpt.*Monod;                                                 % (mol_i/Lr h)
r = r.*(1 + z_r);                                                          % Selecting the reaction for which energetics is being evaluated
R.rm.r = r;
r = r*(pOp.Vr/pOp.Vx);                                                     % (mol_i/Lx h)
R.rm.r = r;

%% Active transport rate determination
my_Act_trans();

%% Energetics determination
indexATP = strcmp(rm.rmNames,'ATPsynthase');
r(indexATP) = my_energetics();
R.rm.r = r;

%%
    function rt = my_Act_trans()
        %% function my_Act_trans
        % rt = my_Act_trans()
        % returns structure R updating:
        % - active transport rates
        % - substrate-related transport rates
        % Inputs:   No inputs because it is a nested function
        % Outputs:
        %           rt           sum of active and passive transport rates
        %% Loading variables
        rmTr = R.rmTr;
        DGtr = R.AlgSt.DGtr;
        numrt = length(DGtr);
        K_act = R.pKt.KTr.K_actV(1:numrt);
        rt_diff  = rmTr.rt_dif(1:numrt);
        StTrans = St.StV(1:numrt);
        M_act = R.pKt.KTr.M_actV(1:numrt);
        z_rt = Opt.z_rt(1:numrt);
        
        %% Active transport rates determination
        
        rt_act = K_act*St.X.*((StTrans)./((M_act)+StTrans));               % (mol_i/Lr h)
        rt_act = rt_act.*(1 + z_rt);
        
        switchFunction = (tanh(R.pOp.Gan*(StTrans-0.01))+1)/2;
        
        reac=(R.rm.stoM(1:numrt,:)*R.rm.r)/(pOp.Vr/pOp.Vx).*sign(K_act);   % (mol_i/Lr h)
        A=reac;
        A=max(A,0);
        
        rt_act = switchFunction.*(A) + (1-switchFunction).*rt_act;
        
        %% Global transport rates determination
        rt  = rt_diff + rt_act;
        rt = [rt;rmTr.rt_dif(numrt+1:end)];
        aux = [zeros(numrt,1);ones((length(rmTr.rt_dif)-numrt),1)];
        rt = (1-aux).*rt*(pOp.Vr/pOp.Vliq) + aux.*rt;                      % (mol_i/Lliq h)
        
        %% Variables update
        R.rmTr.rt_act = rt_act;
        R.rmTr.rt = rt;
        
    end
    function hOpt = my_regulation()
        %% function f_my_regulation
        % hOpt =f_my_regulation(R,t)
        % returns structure R updating:
        % - reaction rates of Glycolysis, Anabolism, Decay and  Na+-pump
        % Inputs:   There are no inputs because it is a nested function
        % Outputs:
        %           hOpt      Rates of Glycolysis, Anabolism, Decay and Na+ pump
        
        hOpt = ones(rm.num_r,1);
        %% Na+ Pump
        index = strcmp(rm.rmNames, 'Na_Pump');
        hOpt(index) = 0;
        
        %% Anabolism
        index = strcmp(rm.rmNames, 'Anab');
        hOpt(index) = 0;
        
        index = strcmp(rm.rmNames, 'AnabProt');
        hOpt(index) = 0;
        
        %% Decay
        index = strcmp(rm.rmNames,'DecayGlu');
        hOpt(index) = 0;
        
        index = strcmp(rm.rmNames,'DecayProt');
        hOpt(index) = 0;
        
        %% Updating variables
        R.rm.hOpt = hOpt;
        R.rm.f = f;
    end
    function rATPsyn = my_energetics()
        %% function my_energetics
        % rATPsyn = my_energetics()
        % returns: pmf-related and maintenance ATP rate
        % Inputs:  There are no inputs because it is a nested function
        %
        % Outputs:
        %          rATPsyn      pmf-related ATP rate
        
        %% Loading variables
        rm = R.rm;
        f = R.rm.f;
        rt_act = R.rmTr.rt_act*(pOp.Vr/pOp.Vx);                            % Active transport rate epxress in mol_i/Lx h
        n_r = rm.num_r;
        c = rm.coupled; c_r=c(1:n_r); c_rt=c(n_r+1:end);                   % Which reactions are related with pmf
        DGr = AlgSt.DGr; DGtr = AlgSt.DGtr; f_r = AlgSt.f.f_r; f_rt = AlgSt.f.f_rt; DGmin = pOp.DGmin; pmf = pOp.pmf;
        [Hup_r, Hdown_r,Hup_rt, Hdown_rt]  = translocations();
        
        pos_ATP=strcmp(rm.rmNames,'ATPsynthase');
        H_atp = Hdown_r(pos_ATP);                                          % pmf-ATP equivalency
        %H_atp = 3;
        rc = [rm.r; rt_act];                                               % (mol_i/Lr h)
        
        aux_r= sign(rm.r) ~= sign(f_r);
        %H_trans_r=min(((1-aux_r).*Hup_r - aux_r.*Hdown_r),1).*c_r;
        H_trans_r=(1-aux_r).*Hup_r - aux_r.*Hdown_r.*c_r;
        
        aux_rt= sign(rt_act) ~= sign(f_rt);
        H_trans_rt=((1-aux_rt).*Hup_rt - aux_rt.*Hdown_rt).*c_rt;
        
        H_trans = [H_trans_r; H_trans_rt];                                 % # of protons exchanged in each reaction (mol H+/reac)
        
        v_atp = sum(H_trans.*abs(rc));                                     % H+ rate (mol H+/Lx h)
        rATPsyn = v_atp/H_atp;                                             % ATP rate (mol ATP/Lr h)
        
        % Contribution of acidic, basic and pmf to the global ATP expenditure (for results analysis only)
        pos_NH4Tr= strcmp(R.rmTr.rmTrNames,'NH3Tr');
        pos_NH4Tr= pos_NH4Tr(1:end-2);
        pos_AcTr = pos_NH4Tr==0;
        R.e.pmf = sum(H_trans(1:n_r).*abs(rc(1:n_r)))./H_atp;
        R.e.Ac = sum(H_trans_rt(pos_AcTr).*abs(rt_act(pos_AcTr)))./H_atp;
        R.e.NH4 = sum(H_trans_rt(pos_NH4Tr).*abs(rt_act(pos_NH4Tr)))./H_atp;
        
        
        %% MAINTENANCE
        m_atp = 0;
        
        rATPsyn = rATPsyn - m_atp;
        
        function [Hup_r, Hdown_r,Hup_rt, Hdown_rt] = translocations()
            
            % Num. of translocated protons for the reactions with DG < 0
            Hup_r=(tanh(10*(DGmin-DGr-c_r.*pmf))+1)/2.*c_r;
            %Hup_r = fix((DGmin-sign(f_r).*DGr)./(pmf.*c_r));
            Hup_r(isinf(Hup_r)) = 0;
            Hup_rt = ((DGmin-sign(f_rt).*DGtr)/pmf);
            
            % Num. of translocated protons for the reactions with DG > 0
            Hdown_r = ceil(((f_r==0).*DGr-DGmin)/pmf);
            Hdown_rt = ((-sign(f_rt).*DGtr-DGmin)/pmf);
        end
    end
end