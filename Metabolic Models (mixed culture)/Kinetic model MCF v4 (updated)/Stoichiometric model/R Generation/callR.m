
clc
clear all

% idR = input('Introduce idR number [1]: ');
% if isempty(idR)
%     idR = 1;
% end

idR=1;
%%%%%%%%%INICIALIZATION...

% Reading the model parameters from Excel file
fprintf('\n> LOADING AND CREATING MODEL STRUCTURE AND PARAMETERS...')
[R, Opt] = loadModelXlsx_PROT(idR);
fprintf('\n... MODEL STRUCTURE AND PARAMETERS LOADED >>>>>\n')

% Executed function to get feed variables
R = my_feeding(0, idR, R);

% At t = 0, efluent composition is the same as influent composition
R(idR).Eff = R(idR).Inf_;
R(idR).Eff.EffV = R(idR).Eff.InfV;
R(idR).Eff = rmfield(R(idR).Eff, 'InfV');

R(idR).St.NOut = R(idR).St.numSt + 2 + 3*R(idR).rm.num_r + 2*R(idR).rmTr.num_tr + 2;

%%% Initial volume of microorganisms
pOp = R(idR).pOp;
St = R(idR).St;
Xt = sum(St.StV(strcmp((St.Phase),'S')));
pOp.Xt = Xt;
Vx = (pOp.Xt*pOp.Vr)/pOp.rho;
pOp.Vx = Vx;
pOp.Vliq = pOp.Vr - Vx;
pOp.dXt_dt = 0;

R(idR).pOp = pOp; 
R(idR).St = St;
R(idR).TolR = 1e-7;
R(idR).TolA = 1e-6;

save R.mat R
save Opt.mat Opt
fprintf('\n...>>>>> R.mat and Opt.mat LOADED READY FOR SIMULATIONS\n')
