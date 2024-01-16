%% Simulate anaerobic digestion reactor based on ADM model as implemented in the benchmark simulation model no. 1 (BSM1)
% Rebeca Gonzalez-Cabaleiro. University of Santiago de Compostela. Spain
% Modified by Alberte Regueira
% Last modification 
%Please contact alberte.regueira@usc.es if you intend to use this code.

tic
clc
clear
close all %Delete figures

load R.mat       %Load the content of R.mat
load Opt.mat

R.pOp.Ana=5;
R.pOp.Gan=5e3;
Opt.electronSource=0;
Opt.eSource=-0.2;
R.pOp.control=[1e-3]; 

x0 = R.St.StV;
vec = R.St.act_states;
x0=x0(vec==1);

nneg = 1:length(x0);      %index vector equal to from 1 to the number of states (step +1). It represents the states that cannot become negative
tspan = 0:1/10:100;
options = odeset('OutputFcn', @(t, x, flag) out(t, x, flag), 'RelTol', R.TolR,'AbsTol', R.TolA,'NonNegative', nneg, 'MaxStep', 0.1); %options for the ode solver
clc
fprintf('\n> MODEL RUNNING >>>>>')

[T,Y] = ode15s(@(t, x) my_model(t, x, R), tspan, x0, options);  %Integrate the system of ode in function my_model, during time defined at tspan, with initial states x0

M = excel(R,Opt);

save M.mat M

