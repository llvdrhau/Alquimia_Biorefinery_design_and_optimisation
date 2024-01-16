function Qgas=f_my_gasPhase(states,compounds, parameters)
%% Returns the characteristics of the gas phase to integrate in the mass balances

% Miguel Mauricio Iglesias. University of Santiago de Compostela. Spain
% November 2014. Please contact miguel.mauricio.iglesias@gmail.com if you
% intend to use this code.

parValues = parameters.reactorPar;
parAbb = parameters.reactorNames;
compoundAbb = compounds.Abb;
Mw = compounds.MW;

Sgas_h2 = states(strcmp(compoundAbb,'Sgas_h2'));    %H2 in gas headspace
Sgas_co2 = states(strcmp(compoundAbb,'Sgas_co2'));  %CO2 in gas headspace
Sgas_ch4 = states(strcmp(compoundAbb,'Sgas_ch4'));  %CH4 in gas headspace
kp = parValues(strcmp(parAbb,'kp'));

R = parValues(strcmp(parAbb,'R'));     %J/K/kmol
T = parValues(strcmp(parAbb,'T'));  %K

Patm = 1.013e5;   % Pa
Ph2o = 0.0313*exp(5290*(1/298-1/T));     %Pa


Ph2 = Sgas_h2*R*T/Mw(strcmp(compoundAbb,'Sh2'));      %Pa
Pch4 = Sgas_ch4*R*T/Mw(strcmp(compoundAbb,'Sch4'));   %Pa
Pco2 = Sgas_co2*R*T/Mw(strcmp(compoundAbb,'Stic'));   %Pa
Pgas = Ph2 + Pch4 + Pco2 + Ph2o;
Qgas = kp*(Pgas - Patm)*Pgas/Patm;