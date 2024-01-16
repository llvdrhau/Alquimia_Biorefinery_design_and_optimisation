function R = my_stoichiometry(R)
%% function my_stoichiometry
% R = my_stoichiometry(R)
% returns: reaction-related rates of each compound
% Inputs:   
%           R           parameter structure
% Outputs:
%           R           parameter structure

%% Loading varaibles
St = R.St;
pOp = R.pOp;
stoM = R.rm.stoM;
r = R.rm.r;
trpM = R.rmTr.trpM;
rt = R.rmTr.rt;

%% Change of base (from reaction base to compound base)

NR = stoM*r; 
aux = ones(length(NR),1) + strcmp(St.Phase, 'S')*((pOp.Vx/pOp.Vr)-1) + strcmp(St.Phase, 'tR')*((pOp.Vx/pOp.Vliq)-1);
NR = aux.*NR;                                                              % Modify units for BM to mol_i/Lr·h, and for extracellular compounds to mol_i/Lliq·h
NR(strcmp(St.StNames,'Ce_Na')) = 0;                                        % Extracelleluar Na+ rate is 0


NRt = trpM*rt;                                                           
aux = ones(length(NRt),1) + strcmp(St.Phase, 'tC')*((pOp.Vliq/pOp.Vx)-1) + strcmp(St.Phase, 'tG')*((pOp.Vliq/pOp.Vgas)-1);
NRt = aux.*NRt;                                                            % Modify units for intracellular-transportable compounds to mol_i/Lx·h and for gas compounds to mol_i/Lgas·h 


%% Variable update
R.NR = NR;
R.NRt = NRt;