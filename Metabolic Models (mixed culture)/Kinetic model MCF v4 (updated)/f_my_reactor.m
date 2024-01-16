function dstatesdt = f_my_reactor(t,states,compounds,parameters, stoiMatrix,pH,Dliq)

% This function calculates the derivative of states with respect to time in
% a continuos stirred tank reactor (CSTR)
% *********************************************************************** %
% INPUTS:
% states for each time
% dilution time
% feeding
% compounds name
% kinetic parameters
% stoichiometric matrix
% pH
% *********************************************************************** %
% OUTPUT:
% derivative of states with respect to time
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


feed=compounds.feed; %feeding concentrations;
Vgas=parameters.reactorPar(strcmp(parameters.reactorNames,'Vgas'));


%% States reaction rate function

states(states<0)=0;          

r= f_my_kinetics(states,compounds,parameters, stoiMatrix, pH); % reaction rate
rt=f_my_mass_transfer(states,compounds,parameters,pH); % transfer rate

flagLiquid = compounds.flagLiquid;
flagGas = compounds.flagGas;
D = states*0;
D(flagLiquid==1) = Dliq;

if any(flagGas)
    [Qgas]=f_my_gasPhase(states,compounds, parameters);
    D(flagGas) = Qgas/Vgas;
end


%% States balance in a CSTR
dstatesdt=D.*feed-D.*states+r+rt;

end

