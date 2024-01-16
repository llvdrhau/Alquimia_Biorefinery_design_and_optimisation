function [ M ] = excel(R,Opt)

x=R.St.StV;

%% my moodel
RUTA=Opt.RUTA;
TRANSP=Opt.TRANSP;
t0=0;
t=0;

St = R.St;                     %Assigns value to states
pOp = R.pOp;                   %Operation parameters
iNae= strcmp(St.StNames,'Ce_Na');
iNai= strcmp(St.StNames,'Ci_Na');
x(iNae) = x(iNai);                    %Imposes the extracellullar and intracellular sodium concentration to be equal
auxCC = (x <= pOp.Clim);            %Is any state less than the concentration limit
x = auxCC.*pOp.Clim + (1-auxCC).*x; % The state is equal to the minimum concentration if the concentration falls below
% Update of the state variables
StV = x; %StV= Values of states.

aux = St.StNames;
for i=1:(St.numSt)
    St.(char(aux(i))) = x(i);       %Assigns the corrected value of each state to the corresponding field in the structure St. TO BE CONSIDERED OR RECODED
end

StVchr = [StV(1:St.pos_CitC)                    % Intracellular states that vary (1-24 in StV == x)
    1                                   % This 1 is equal to the internal water concentration. It is used linked with matrix Keq and chrM in pH functions
    StV((St.pos_CitC+1):St.pos_Ci)      % Intracellular moieties (25-35 in StV)
    StV(St.pos_S:St.pos_Ce)             % Biomass and extra cellular quantities (36-60 in StV)
    1                                   % This 1 is equal to the external water concentration. It is used linked with matrix Keq and chr in pH functions.
    StV((St.pos_Ce+1):end)];            % Gas states (61-63 in StV)

% Update in the global variable the states for use by other functions.
R.St = St;
R.St.StV = StV;
R.St.StVchr = StVchr;

 R = my_algebraics(R,t);
    
    % Microorganisms volume
    Xt = sum(R.St.StV(strcmp((R.St.Phase),'S'))); %Biomass concentration
    pOp.Xt = Xt;
    Vx = (pOp.Xt*pOp.Vr)/pOp.rho; %Biomass volume (L)
    pOp.Vx = Vx;
    % Liquid volume
    pOp.Vliq = pOp.Vr - Vx;
    
       R.pOp = pOp; %Update of pOp in R.
    % Calculation of the transport rates
    R = my_transport(R);

%% my_z
R.flagOpt = 1; %CHANGE OF FLAG
RUTAtot = RUTA;
TRANSPtot = TRANSP;
index_nz=Opt.index_nz;
pos_nz=Opt.pos_nz;
x0 = St.StV;
act_states = R.St.act_states;
x0=x0(act_states==1);

dATP = zeros(1, size(RUTAtot,1)); 
dNADH = zeros(1, size(RUTAtot,1)); 
dCO2 = zeros(1, size(RUTAtot,1)); 
M.e = zeros(size(RUTAtot,1),3);
R.Opt = Opt;
RR = struct(R);
nz = size(RUTAtot,1); %24. Number of pyruvate-related reactions that happen.

for i = 1: nz
    RR(i) = R;
    RUTA = RUTAtot(i,:)'; %i row of RUTA (i.e reaction number of the i of pyruvate-consuming reactions)
    TRANSP = TRANSPtot(i,:)'; %First row of TRANSP (i.e transp reactions used in the number i pyruvate-consuming reaction).
    RR(i).Opt.z_r=RUTA-1; %0=Happening reaction -1=Either non-happening or deactivated reactions.
    RR(i).Opt.z_rt=TRANSP-1; %0=Happening reaction -1=Either non-happening or deactivated reactions.
end
%The difference among the different RR is just z_r and z_rt. Therefore the differences in line 42 will be in the execution of my_kinetics.
f_r=RUTAtot*R.AlgSt.f.f_r;

pos1 = find(strcmp(St.StNames,'Ce_Glu'));
posend = find(strcmp(St.StNames,'Ce_Hist'));
feeding = St.StV([pos1 pos1+2:posend]);
no_feeding = find(feeding<1e-8);
% no_feeding(length(no_feeding)+1) = 3;
desact = ones(length(dATP),1);
for i=1:length(no_feeding)
    if no_feeding(i)==1
        f_r(1:index_nz(1)) = 0;
        desact(1:index_nz(1)) = 0;
    else
        f_r(index_nz(no_feeding(i)-1)+1:index_nz(no_feeding(i))) = 0;
        desact(index_nz(no_feeding(i)-1)+1:index_nz(no_feeding(i))) = 0;
    end
end
optimise = find(desact==1);

for k = 1:length(optimise)
    j=optimise(k);
    [~, RS] = my_model(t0, x0, RR(j));     
    
    if Opt.electronSource
        eSource = Opt.eSource;
    else
        eSource=0;
    end
    
    a = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_NADH'), 4:end).*RS.rm.r(4:end)');  %r_NADH
    b = -(sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_Pyr'), 4:end).*RS.rm.r(4:end)'));%-r_Pyr
    
    dATP(j) = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_ATP'), 4:end).*RS.rm.r(4:end)')+b;
    dCO2(j) = sum(RS.rm.stoM(strcmp(RS.St.StNames, 'Ci_CO2'), 4:end).*RS.rm.r(4:end)');
    dNADH(j) = (1+eSource)*b + a;
    
    M.e(j,1)=RS.e.pmf;
    M.e(j,2)=RS.e.Ac;
    M.e(j,3)=RS.e.NH4;
        
end


dATP((dATP<-100))=-100;
dNADH(isinf(dNADH)) = -100;
dATP(isnan(dATP)) = 0;
dNADH(isnan(dNADH)) = 0;
dCO2(isnan(dCO2)) = 0;

%f_r=RUTAtot*R.AlgSt.f.f_r;

f_r(index_nz)=1;

[z,EXITFLAG,LAMBDA,deltaZ] = f_quadraticProgram(Opt.z,dATP,dNADH,nz,index_nz,f_r);

zmax=max(deltaZ);
fprintf('Max value in deltaZ is: %f \n',zmax)

%% my_kinetics

Kr = R.pKt.Kin.KrV;    z_r = Opt.z_r;
St = R.St;             pOp = R.pOp;
rm = R.rm;             AlgSt = R.AlgSt;
stoM = full(R.rm.stoM); 
K = R.pKt.Kin.K;

if Opt.electronSource
    eSource=Opt.eSource;
else
    eSource=0;
end

f = AlgSt.f.f_r; 
f = f*St.X;                                                                % (mol X/Lr)
R.rm.f = f;

%% 
R = f_my_regulation(R,Opt,t);                                                  % Determines rate of Anabolism, Decay and Na+-pump
h=R.rm.h;
f=R.rm.f;

%% Loading Monod terms
pos_Glu=Opt.pos_Glu;
% pos_Arg=Opt.pos_Arg;
% pos_Ala=Opt.pos_Ala;
% pos_Asp=Opt.pos_Asp;
% pos_Lys=Opt.pos_Lys;
% pos_Glut=Opt.pos_Glut;
% pos_Ser=Opt.pos_Ser;
% pos_Thr=Opt.pos_Thr;
% pos_Cys=Opt.pos_Cys;
% pos_Gly=Opt.pos_Gly;
% pos_Prol=Opt.pos_Prol;
% pos_Vali=Opt.pos_Vali;
% pos_IsoL=Opt.pos_IsoL;
% pos_Leu=Opt.pos_Leu;
% pos_Meth=Opt.pos_Meth;
% pos_GluM=Opt.pos_GluM;
% pos_AspG=Opt.pos_AspG;
% pos_Hist=Opt.pos_Hist;
Monod = Opt.Monod;

%Monod(pos_Glu)=(St.Ce_Glu/(1e-3+St.Ce_Glu));
% Monod=ones(size(Kr))*Monod_Glu;
% Monod(1:pos_EMP)=1;        % These reactions have already been determined: Anabolism, Decay and Glycolysis. These reactions are not limited by auxGlu
% Monod(pos_Arg)=Opt.Monod_Arg;
% Monod(pos_Ala)=Opt.Monod_Ala;
% Monod(pos_Asp)=Opt.Monod_Asp;
% Monod(pos_Lys)=Opt.Monod_Lys;
% Monod(pos_Glut)=Opt.Monod_Glut;
% Monod(pos_Ser)=Opt.Monod_Ser;
% Monod(pos_Thr)=Opt.Monod_Thr;
% Monod(pos_Cys)=Opt.Monod_Cys;
% Monod(pos_Gly)=Opt.Monod_Gly;
% Monod(pos_Prol)=Opt.Monod_Prol;
% Monod(pos_Vali)=Opt.Monod_Vali;
% Monod(pos_IsoL)=Opt.Monod_IsoL;
% Monod(pos_Leu)=Opt.Monod_Leu;
% Monod(pos_Meth)=Opt.Monod_Meth;
% Monod(pos_GluM)=Opt.Monod_GluM;
% Monod(pos_AspG)=Opt.Monod_AspG;
% Monod(pos_Hist)=Opt.Monod_Hist;
% Monod(find(strcmp(rm.rmNames,'NADH > H2')):end) = 1;  %Moiety related reactions (index 34:39) and Na Pump (index 40)
%% Reaction rate determination
r = f.*K.*Kr.*h.*Monod;  
r = r*(pOp.Vr/pOp.Vx);
r = Opt.RUTA*r;

pos_ATP=strcmp(St.StNames,'Ci_ATP');
matriz=(Opt.RUTA*stoM')';
rATP=(matriz(pos_ATP,:).*r')';
coeff_ATP=matriz(pos_ATP,:)';

%% Outputs
M.dATP=dATP';
M.dNADH=dNADH';
M.rATP=rATP;
M.coeff_ATP=coeff_ATP;
M.z=z;
M.deltaZ=deltaZ;
M.f_r=f_r;
M.LAMBDA=LAMBDA;
M.r=r;

output(:,1)=1:length(z);
output(:,2)=f_r;
output(:,3)=dNADH';
output(:,4)=dATP';
output(:,5)=rATP;
output(:,7:9)=M.e;
output(:,6)=output(:,8)+output(:,9);
output(:,10)=z;
output(:,11)=coeff_ATP;

M.output=output;
M.R=R;
M.Opt=Opt;

save M.mat M

%% ATP

ATP_yield{1,1}='ATP/Lx·h';
ATP_yield{1,2}='ATP/h';
ATP_yield{1,3}='BM g/L';

ATP_yield{2,1}=sum(stoM(pos_ATP,6:end).*R.rm.r(6:end)');
ATP_yield{2,2}=ATP_yield{2,1}*R.pOp.Vx;
ATP_yield{2,3}=R.pOp.Xt*24.6;

M.ATP_yild=ATP_yield;
end

