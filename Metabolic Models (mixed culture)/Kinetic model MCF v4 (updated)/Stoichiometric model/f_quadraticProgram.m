function [newZ,EXITFLAG,LAMBDA,deltaZ] = f_quadraticProgram(oldZ,dATP,dNADH,nz,index_nz,f_r)
%% function f_my_quadraticProgram
% [newZ,EXITFLAG,LAMBDA,deltaZ] = f_quadraticProgram(oldZ,dATP,dNADH,dCO2,nz,index_nz,f_r)
% returns the reaction optimisation vector and monitoring optimisation parameters
% Inputs:   oldZ         Reaction optimised vector of the previous time step
%           dATP         ATP rate of the reactions
%           dNADH        NADH rate of the reactions
%           dCO2         CO2 rate of the reactions
%           nz           Number of eligible reactions
%           index_nz     Position of the first reaction of the different substrates
%           f_r          Feasability vector
% Outputs:
%           newZ         New reaction optimised vector
%           EXITFLAG     Flag of the optimiser
%           LAMBDA       Set of Lagrangian multipliers
%           deltaZ       Variation of reaction optimised vector


fMax = -dATP; 
weightMoves = 3;                                                           % Weight factor 
H = eye(nz)*weightMoves;    
f = fMax';                
%% Inequality constraint
A=[];                    
b = [];                                                               % Inequality constraint (CO2 cannot be consumed)

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
beq = [-dNADH*oldZ;zeros(length(index_nz),1);-oldZ(x)];   

% Aeq(1,:)=[];
% beq(1)=[];
%% Boundaries
lowerBoundary = -oldZ;       
upperBoundary = (1-oldZ);   
op = optimset('Display','off','Algorithm','interior-point-convex');

%% Optimiser
[deltaZ,~,EXITFLAG,~,LAMBDA] = quadprog(H,f,A,b,Aeq,beq,lowerBoundary,upperBoundary,[],op);
newZ = oldZ+deltaZ;  
end

