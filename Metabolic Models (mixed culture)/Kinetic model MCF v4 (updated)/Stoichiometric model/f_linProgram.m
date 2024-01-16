function [z,EXITFLAG,LAMBDA]=f_linProgram(oldZ,dATP,dNADH,index_nz,f_r)
%% function f_linProgram
% [z,EXITFLAG]=f_linProgram(dATP,dNADH,dCO2,f_r,oldZ,index_nz)
% returns the reaction optimisation vector and monitoring optimisation parameters
% Inputs:   oldZ         Reaction optimised vector of the previous time step
%           dATP         ATP rate of the reactions
%           dNADH        NADH rate of the reactions
%           dCO2         CO2 rate of the reactions
%           index_nz     Position of the first reaction of the different substrates
%           f_r          Feasability vector
% Outputs:
%           z            New reaction optimised vector
%           EXITFLAG     Flag of the optimiser
%           LAMBDA       Set of Lagrangian multipliers


fMax = -dATP;  

%% Inequality constraint
A = []; b = [];                                                          % Inequality constraint (CO2 cannot be consumed)

%% Equality constraints
Aeq(1,:) = dNADH;                                                          % NADH conservation
Aeq(2,1:index_nz(1)) = 1;                                                  % Reaction optimised vector must sum 1 for each substrate
for  k =1:length(index_nz)-1
   Aeq(k+2,(index_nz(k)+1):index_nz(k+1))=1;
end

x=find(f_r<1e-3);                                                          % Endergonic reactions are not possible
[ny,~]=size(Aeq);
for k=1:length(x)
    Aeq(ny+k,x(k))=1;
end
beq=[0;ones(length(index_nz),1);zeros(length(x),1)];

%% Boundaries
lowerBoundary = zeros(size(dATP));
upperBoundary = ones(size(dATP));
initialGuess=oldZ;


%% Optimiser
op = optimset('Display','off');
[z,~,EXITFLAG,~,LAMBDA]=linprog(fMax,A,b,Aeq,beq,lowerBoundary,upperBoundary,initialGuess,op);
end