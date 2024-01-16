function rt=f_my_mass_transfer(states,compounds,parameters,pH)
%% Gas-liquid transfer between the liquid hold-up and the headspace
% Miguel Mauricio Iglesias. University of Santiago de Compostela. Spain
% November 2014. Please contact miguel.mauricio.iglesias@gmail.com if you
% intend to use this code.

rt = states*0;

parValues = parameters.reactorPar;
parAbb = parameters.reactorNames;
compoundAbb = compounds.Abb;
Mw =compounds.MW;

H = 10^-pH;
rt = states*0;

%% Identification of compounds and parameters
iTic = strcmp(compoundAbb,'Stic');
iH2 = strcmp(compoundAbb,'Sh2');
iCH4 = strcmp(compoundAbb,'Sch4');

Stic = states(iTic);      %TIC
Sh2 = states(iH2);        %H2
Sch4 = states(iCH4);      %CH4
Sgas_h2 = states(strcmp(compoundAbb,'Sgas_h2'));    %H2 in gas headspace
Sgas_co2 = states(strcmp(compoundAbb,'Sgas_co2'));  %CO2 in gas headspace
Sgas_ch4 = states(strcmp(compoundAbb,'Sgas_ch4'));  %CH4 in gas headspace

Kac_tic = parValues(strcmp(parAbb,'Kac_tic'));
kLa = parValues(strcmp(parAbb,'kLa'));
H_co2 = parValues(strcmp(parAbb,'H_co2'));
H_ch4 = parValues(strcmp(parAbb,'H_ch4'));
H_h2 = parValues(strcmp(parAbb,'H_h2'));

R = parValues(strcmp(parAbb,'R'));     %J/K/kmol
T = parValues(strcmp(parAbb,'T')); %K

%% Stripping of CO2, H2 and CH4

Vliq= parValues(strcmp(parAbb,'Vliq'));
Vgas = parValues(strcmp(parAbb,'Vgas'));

% Partial pressure (Pa)
Ph2 = Sgas_h2*R*T/Mw(iH2);      %Pa
Pch4 = Sgas_ch4*R*T/Mw(iCH4);   %Pa
Pco2 = Sgas_co2*R*T/Mw(iTic);   %Pa

Sco2 = Stic*(H/(H+Kac_tic));

% Transfer
rt(iH2) = kLa*(Mw(iH2)*H_h2*Ph2-Sh2);         %kgCOD/d/m3
rt(iCH4) = kLa*(Mw(iCH4)*H_ch4*Pch4-Sch4);    %kgCOD/d/m3
rt(iTic) = kLa*(Mw(iTic)*H_co2*Pco2-Sco2);    %kgCOD/d/m3

rt(strcmp(compoundAbb,'Sgas_h2')) = -rt(iH2)*Vliq/Vgas;       %kgCOD/d/m3
rt(strcmp(compoundAbb,'Sgas_ch4')) = -rt(iCH4)*Vliq/Vgas;     %kgCOD/d/m3
rt(strcmp(compoundAbb,'Sgas_co2')) = -rt(iTic)*Vliq/Vgas;     %kgCOD/d/m3

end
