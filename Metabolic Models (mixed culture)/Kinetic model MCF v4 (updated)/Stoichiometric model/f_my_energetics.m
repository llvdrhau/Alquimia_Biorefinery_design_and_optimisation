function rATPsyn = f_my_energetics(R)
%% function f_my_energetics
% rATPsyn = f_my_energetics(R)
% returns: pmf-related and maintenance ATP rate
% Inputs:   
%           R           parameter structure
% Outputs:
%          rATPsyn      pmf-related ATP rate

%% Loading varaibles
AlgSt = R.AlgSt;
pOp = R.pOp;
St = R.St;  
rm = R.rm;
rt_act = R.rmTr.rt_act*(pOp.Vr/pOp.Vx);                                    % Active transport rate epxress in (mol_i/Lx h)
n_r = rm.num_r;
c = rm.coupled; c_r=c(1:n_r); c_rt=c(n_r+1:end);                           % Which reactions are related with pmf
DGr = AlgSt.DGr; DGtr = AlgSt.DGtr; f_r = AlgSt.f.f_r; f_rt = AlgSt.f.f_rt; DGmin = pOp.DGmin; pmf = pOp.pmf;

%% pmf-related ATP rate
[Hup_r, Hdown_r,Hup_rt, Hdown_rt] = translocations();                      

pos_ATP=strcmp(rm.rmNames,'ATPsynthase');
H_atp = Hdown_r(pos_ATP);                                                  % pmf-ATP equivalency
%H_atp = 3;
rc = [rm.r; rt_act];                                                       % (mol_i/Lx h)

aux_r= sign(rm.r) ~= sign(f_r); 
%H_trans_r=min(((1-aux_r).*Hup_r - aux_r.*Hdown_r),1).*c_r;
H_trans_r=(1-aux_r).*Hup_r - aux_r.*Hdown_r.*c_r;

aux_rt= sign(rt_act) ~= sign(f_rt);
H_trans_rt=((1-aux_rt).*Hup_rt - aux_rt.*Hdown_rt).*c_rt;

H_trans = [H_trans_r; H_trans_rt];                                         % # of protons exchanged in each reaction (mol H+/reac)

v_atp = sum(H_trans.*abs(rc));                                             % H+ rate (mol H+/Lx h)
rATPsyn = (v_atp/H_atp);                                                   % ATP rate (mol ATP/Lx h)

%% Maintenance-reltaed ATP rate
DGrATP = AlgSt.DGr(strcmp(rm.rmNames,'ATPsynthase'));
m_atp = (4.5/DGrATP)*St.X;
rATPsyn = rATPsyn - m_atp;                                                 % ATP rate (mol ATP/Lx h)

    function [Hup_r, Hdown_r,Hup_rt, Hdown_rt] = translocations()
        
        % Num. of translocated protons for the reactions with DG < 0
        %Hup_r = fix((DGmin-sign(f_r).*DGr)./(pmf.*c_r));
        Hup_r=(tanh(10*(DGmin-DGr-c_r.*pmf))+1)/2.*c_r;
        Hup_r(isinf(Hup_r)) = 0;
        Hup_rt = ((DGmin-sign(f_rt).*DGtr)/pmf);
        
        % Num. of translocated protons for the reactions with DG > 0
        Hdown_r = ceil(((f_r==0).*DGr-DGmin)/pmf);
        Hdown_rt = ((-sign(f_rt).*DGtr-DGmin)/pmf);
    end
end