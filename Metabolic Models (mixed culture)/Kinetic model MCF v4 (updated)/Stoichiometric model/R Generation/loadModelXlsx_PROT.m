% loadModelXls.m -Load all the parameters from Excel sheet
function [R, Opt] = loadModelXlsx_PROT(idR)
warning off all

% If only one reactor is used is set to '' and no extra for the names of the xls sheets is needed
if idR==1,  id = '';    else    id = num2str(idR);   end
route = char('Excel\myModelStructure.xlsx');

% GENERAL MODEL PARAMETERS from the Excel file
fprintf('\n> Loading and creating GENERAL MODEL PARAMETERS...')

[OperatParam, pOpNames]   = xlsread(route, strcat('OperatParam ', id));
aux = '';
for i=1:length(OperatParam)
    % Building the string command for the structure, last without ', '
    if i<length(OperatParam),
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)),', ');
    else
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)));
    end
end
aux = char(aux);
eval(strcat('pOp = struct(', aux, ');'));

% Structure with the THERMODYNAMIC PARAMETERS
[ThParam, pThNames]  = xlsread(route, strcat('ThermoParam', id));
ThParam = ThParam(:,1:5);

aux = '';
for i=1:size(ThParam,1)
    % Building the string command for the structure, last without ', '
    if i<size(ThParam,1)
        aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:), ');
    else
        aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:)');
    end
end
aux = char(aux);
eval(strcat('pTh = struct(', aux, ');'));
pTh.pThNames = pThNames(:,1);
pTh.Type = pThNames(:,size(pThNames,2));

idnSp = prod(1*(char(pTh.Type) == repmat('Met'|'Aux'|'Sub'|'FMe'|'Bio'|'Gas',size(pThNames,1),1)),2) == 1;

idG0f = sum([(prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('G0l',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('G0s',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('G0g',size(ThParam,1),1)),2) == 1)],2);
idH0f = sum([(prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('H0l',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('H0s',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('H0g',size(ThParam,1),1)),2) == 1)],2);
idchr = sum([(prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('Sub',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('Met',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('FMe',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('Aux',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('Bio',size(ThParam,1),1)),2) == 1) ...
    (prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('Gas',size(ThParam,1),1)),2) == 1)],2);
pTh.spcNames= [];
pTh.G0fM = zeros(sum(idG0f),size(ThParam,2));   jG = 0;
pTh.H0fM = zeros(sum(idH0f),size(ThParam,2));   jH = 0;
pTh.chrM = zeros(sum(idchr),size(ThParam,2));   jch = 0;
pTh.spcR = [];
for i=1:size(pThNames,1)
    if  idnSp(i) == 1
        pTh.spcNames = [pTh.spcNames; pThNames(i,2:size(pThNames,2)-1)];
    end
end

for i=1:size(ThParam,1)
    if idG0f(i) == 1
        jG = jG+1;
        pTh.G0fM(jG,:) = ThParam(i,:);
    elseif idH0f(i) == 1
        jH = jH+1;
        pTh.H0fM(jH,:) = ThParam(i,:);
    elseif idchr(i) == 1
        jch = jch+1;
        pTh.chrM(jch,:) = ThParam(i,:);
    end
end
nan_locations = isnan(pTh.G0fM);
pTh.G0fM(nan_locations) = 1e10;
nan_locations = isnan(pTh.chrM);
pTh.chrM(nan_locations) = 0;

% Temperature modified Gibbs free energies of formation
%%%%%%%%%%%%%%%%%%%%pTh.GfM = pTh.H0fM + (pOp.T/298.15)*(pTh.G0fM -
%%%%%%%%%%%%%%%%%%%%pTh.H0fM);!!!!!!!!!!!!!OJOOOOOOOOOO
pTh.GfM = pTh.G0fM;

% Calculation of the acid constants from the Gibbs energy formation matrix provided in thM for the liquid states. The thM matrix provide Gf of
% [Uncharged non hydrated species, (hydrated if aplicable) fully protonated species,
% 1st deprotonated species, 2nd deprotonated species and 3rd deprotonated species]
% K are calculated from exp{(Gºf_prot - Gºf_deprot)/(R·T)}
% Submatrix containing only the Gf of the liquid state variables

idGf_L = prod(1*(char(pTh.Type(1:size(ThParam,1))) == repmat('G0l',size(ThParam,1),1)),2) == 1;
GfM_L = zeros(sum(idGf_L),size(ThParam,2));jGl = 0;
for i=1:size(ThParam,1)
    if idGf_L(i) == 1
        jGl = jGl + 1;
        GfM_L(jGl,:) = pTh.GfM(i,:);
    end
end
Keq = zeros(size(GfM_L,1),4);
Keq(:,1) = exp((pTh.Gf_SH2O(2) + GfM_L(:,1) - GfM_L(:,2))/(-pOp.Rth*pOp.T));
for i=2:(size(GfM_L,2) - 1)
    Keq(:,i) = exp((GfM_L(:,i+1)+(i*1e10*(GfM_L(:,i+1) == 1e10))-(GfM_L(:,i)))/(-pOp.Rth*pOp.T));
end
pTh.GfM_L = GfM_L;
pTh.Keq = Keq;

%% STOICHIOMETRY/REACTION MATRIX and TRANSPORT MATRIX structures
[fullMatrix, rmNamesFull] = xlsread(route, strcat('ReactionMatrixFull', id));
coupled = fullMatrix(end,:)';
fullMatrix = fullMatrix((1:end-1),:);
fullMatrix = sparse(fullMatrix);
StNamesFull = rmNamesFull(2:(end-1),1);
rmNamesFull = rmNamesFull(1,2:end);

aux_r = 1 - (strcmp('Ci_H2O', StNamesFull) + strcmp('Ci_H', StNamesFull) + strcmp('Ce_H2O', StNamesFull) + strcmp('Ce_H', StNamesFull));
aux_c = find(strcmp('GluTr',rmNamesFull) + strncmp(rmNamesFull, 'G_',2));

% Modification of the coupled vector -->
%rm.coupled = coupled(1:(aux_c(2)-1)) + strcmp(rmNamesFull(1:(aux_c(2)-1)),'ATPsynthase')';
rm.coupled = coupled(1:(aux_c(2)-1));
rm.stoM = fullMatrix(aux_r ~= 0, 1:(aux_c(1)-1));
rm.stoMfull = fullMatrix(:, 1:(aux_c(1)-1));
rm.rmNames = rmNamesFull(1:(aux_c(1)-1));
rm.num_r = length(rm.rmNames);


rmTr.trpM = fullMatrix(aux_r~= 0, aux_c(1):end);
rmTr.trpMfull = fullMatrix(:, aux_c(1):(aux_c(2)-1)); % Excluded the liquid-gas transfer reactions (Gibbs calculations)
rmTr.rmTrNames = rmNamesFull(aux_c:end);
rmTr.num_tr = length(rmTr.rmTrNames);

% Inicialization of some variables
rm.h = ones(rm.num_r,1);
rm.h(strcmp(rm.rmNames, 'Na_Pump')) = 0;


%% STATE VARIABLES Structure and initial values
[States, StNames] = xlsread(route, strcat('States',id));

aux = '';
for i=1:length(States)
    % Building the string command for the structure, last without ', '
    if i<length(States)
        aux = strcat(aux, char(39), StNames(i,1), char(39), ', ',num2str(States(i)),', ');
    else
        aux = strcat(aux, char(39), StNames(i,1), char(39), ', ',num2str(States(i)));
    end
end
aux = char(aux);
eval(strcat('St = struct(', aux, ');'));
St.Ci_H2O = 1;
St.Ce_H2O = 1;
St.Ci_H = 1e-7;
St.Ce_H = 10^(-pOp.pH);
St.StNamesFull = StNamesFull;
St.StNames = StNames(:,1);
St.Phase = StNames(:,4);
spcR = States(:,4);
spcN= States(:,5);
St.StV = States(:,1);
St.act_states = States(:,6);

[~,SubsNames] = xlsread(route, strcat('Substrates',id));
St.SubsNames = SubsNames(:,1);
St.nSubs = length(St.SubsNames);
% Number of states
St.numSt = length(St.StV);

for i=1:St.numSt
    if  strcmp(St.Phase(i),'tC') && strcmp(St.Phase(i+1),'mC')
        St.pos_CitC = i;
    elseif  strcmp(St.Phase(i),'mC') && strcmp(St.Phase(i+1),'S')
        St.pos_Ci = i;
        St.pos_Cichr = i+1;
    elseif  strcmp(St.Phase(i),'S') && strcmp(St.Phase(i+1),'tR')
        St.pos_S = i;
        St.pos_Schr = i+1;
    elseif strcmp(St.Phase(i),'tR') && strcmp(St.Phase(i+1),'tG')
        St.pos_Ce = i;
        St.pos_Cechr = i+2;
    end
end
ind_H2O = [find(strcmp('Ci_H2O',StNamesFull));find(strcmp('Ce_H2O',StNamesFull))];
ind_H = [find(strcmp('Ci_H',StNamesFull));find(strcmp('Ce_H',StNamesFull))];
ind_H2O = [(ind_H2O(1) - (ind_H(1)<ind_H2O(1))); (ind_H2O(2) - 1 -(ind_H(2)<ind_H2O(2)))];

St.StVchr = [St.StV(1:(ind_H2O(1)-1));St.Ci_H2O;St.StV(ind_H2O(1):(ind_H2O(2)-2));St.Ce_H2O;St.StV((ind_H2O(2)-1):end)];
St.StchrNames = [St.StNames(1:(ind_H2O(1)-1));'Water';St.StNames(ind_H2O(1):(ind_H2O(2)-2));'Water';St.StNames((ind_H2O(2)-1):end)];
pTh.spcR = spcR;
pTh.spcRfull = [spcR(1:(ind_H2O(1)-1));2;spcR(ind_H2O(1):(ind_H2O(2)-2));2;spcR((ind_H2O(2)-1):end)];
pTh.spcRfull = [pTh.spcRfull(1:(ind_H(1)-1));2;pTh.spcRfull(ind_H(1):(ind_H(2)-2));2;pTh.spcRfull((ind_H(2)-1):end)];% Included the H2O and the H+

pTh.spcN = spcN;

aux = find(strcmp(pTh.Type,'Sub'));
pTh.Path = [pTh.Type((aux):(aux+ind_H(1)-2));'Aux';pTh.Type((aux+ind_H(1)-1):(aux+St.pos_Cichr));pTh.Type(aux:(aux+ind_H2O(1)-1));'Aux';pTh.Type((aux+St.pos_Cichr+1):end)];
%Path = StNames(:,7);
%pTh.Path2 = [Path(1:ind_H(1)-2);'Aux';Path(ind_H(1)-2:St.pos_Ce);'Aux';'Aux';Path(St.pos_Ce+1:end)];
% Calculation of Henry's constants from Gºf values of the soluble gases
vaux = find(strcmp(St.Phase,'tG'));
KhV = zeros(length(St.StV),1);
auxspcR = pTh.spcR;
auxspcR(strcmp(St.StNames,'Ce_CO2')) = 1;
for i=1:length(vaux)
    name = char(St.StNames(vaux(i)));
    name = name(3:end);
    w = strcmp(strcat('Ce_',name),St.StNames);
    Kh.(name) = exp((pTh.(strcat('Gf_S',name))(auxspcR(w)) - pTh.(strcat('Gf_G',name))(auxspcR(vaux(i))))  / (-pOp.Rth*pOp.T));% M/atm
    KhV(vaux(i)) = Kh.(name);
end

pTh.Kh = Kh;
pTh.Kh.KhV = KhV;


%% Structure with the KINETIC PARAMETERS
[KinetParam,  pKtNames] = xlsread(route, strcat('KinetParam',  id));

aux = strcmp('Diffusivity',pKtNames(1,:));
K_dif = KinetParam((isfinite(KinetParam(:,aux))),aux);
K_difNames = pKtNames(2:end,aux);

aux = strcmp('Parameter',pKtNames(1,:));
K = KinetParam((isfinite(KinetParam(:,aux))),aux);
Kin.K = K;

aux = strcmp('Kinetic',pKtNames(1,:));
Kr = KinetParam((isfinite(KinetParam(:,aux))),aux);
KrNames = pKtNames(2:end,aux);

aux = strcmp('Active_Tr',pKtNames(1,:));
K_act = KinetParam((isfinite(KinetParam(:,aux))),aux);
K_actNames = pKtNames(2:end,aux);

aux = strcmp('MonodActive_Tr',pKtNames(1,:));
M_act = KinetParam((isfinite(KinetParam(:,aux))),aux);
M_actNames = pKtNames(2:end,aux);

aux = strcmp('Monod',pKtNames(1,:));
Monod_sat = KinetParam((isfinite(KinetParam(:,aux))),aux);
Monod_Names = pKtNames(2:length(Monod_sat)+1,aux);

aux = '';
for i=1:length(K_dif),
    % Building the string command for the structure, last without ', '
    if i<length(K_dif),
        aux = strcat(aux, char(39), K_difNames(i), char(39), ', ',num2str(K_dif(i)),', ');
    else
        aux = strcat(aux, char(39), K_difNames(i), char(39), ', ',num2str(K_dif(i)));
    end
end
KTr.K_difNames = K_difNames;
KTr.K_difV = K_dif;

aux = '';
for i=1:length(Kr),
    % Building the string command for the structure, last without ', '
    if i<length(Kr),
        aux = strcat(aux, char(39), KrNames(i), char(39), ', ',num2str(Kr(i)),', ');
    else
        aux = strcat(aux, char(39), KrNames(i), char(39), ', ',num2str(Kr(i)));
    end
end
Kin.KrNames = KrNames;
Kin.KrV = Kr;

aux = '';
for i=1:length(K_act),
    % Building the string command for the structure, last without ', '
    if i<length(K_act),
        aux = strcat(aux, char(39), K_actNames(i), char(39), ', ',num2str(K_act(i)),', ');
    else
        aux = strcat(aux, char(39), K_actNames(i), char(39), ', ',num2str(K_act(i)));
    end
end
KTr.K_actNames = K_actNames;
KTr.K_actV = K_act;

for i=1:length(M_act),
    % Building the string command for the structure, last without ', '
    if i<length(M_act),
        aux = strcat(aux, char(39), M_actNames(i), char(39), ', ',num2str(M_act(i)),', ');
    else
        aux = strcat(aux, char(39), M_actNames(i), char(39), ', ',num2str(M_act(i)));
    end
end
KTr.M_actNames = M_actNames;
KTr.M_actV = M_act;

% aux = '';
% for i=1:length(K_),
%     % Building the string command for the structure, last without ', '
%     if i<length(K_),
%         aux = strcat(aux, char(39), K_Names(i), char(39), ', ',num2str(K_(i)),', ');
%     else
%         aux = strcat(aux, char(39), K_Names(i), char(39), ', ',num2str(K_(i)));
%     end
% end
% aux = char(aux);
% eval(strcat('K_ = struct(', aux, ');'));
pKt.Monod_sat = Monod_sat;
pKt.Monod_Names = Monod_Names;
pKt.Kin = Kin;
pKt.KTr = KTr;

%% FEEDING variables and values
[FeedParam, FdNames] = xlsread(route, strcat('FeedProgram',  id));
FdNames = FdNames(2:end);    % Removing "Time" name from the names vector
Time = FeedParam(:,1);
FeedParam = FeedParam(:,2:end);    % Removing "Time" name from the values vector
FdV = FeedParam;
aux = '';
for i=1:length(FdNames),
    % Building the string command for the structure, last without ', '
    if i<length(FdNames),
        aux = strcat(aux, char(39), FdNames(i), char(39), ', FeedParam(:,', num2str(i), '), ');
    else
        aux = strcat(aux, char(39), FdNames(i), char(39), ', FeedParam(:,', num2str(i), ')');
    end
end
aux = char(aux);
eval(strcat('mFd = struct(', aux, ');'));

mFd.FdNames = FdNames;
mFd.Time = Time;
mFd.FdV = FdV;
mFd.Now = 1;
% Influent values at t = 0
InfV = mFd.FdV(1,:);
Inf_.Qgas = InfV(1);
Inf_.Q = InfV(2);
Inf_.InfV = InfV';

for i=1:length(InfV),
    Inf_.(char(mFd.FdNames(i))) = InfV(i);
end

%% OPTIMIZATION Values
[OptimizParam, OptimizNames]   = xlsread(route, strcat('OptimizParam ', id));
aux = '';
for i=1:length(OptimizParam)
    % Building the string command for the structure, last without ', '
    if i<length(OptimizParam),
        aux = strcat(aux, char(39), OptimizNames(i), char(39), ', ',num2str(OptimizParam(i)),', ');
    else
        aux = strcat(aux, char(39), OptimizNames(i), char(39), ', ',num2str(OptimizParam(i)));
    end
end
aux = char(aux);
eval(strcat('Opt = struct(', aux, ');'));

%% Creation of RUTA and TRANSP

Path = pTh.Path(1:(St.pos_Schr));
pos_EMP = find(strcmp(KrNames,'K_EMP'));
sto_matrix = rm.stoMfull(:, pos_EMP:end); %Erase two first reactions: Anab and Decay because they always take place and are calculated elsewhere
num_reactions = rm.num_r-(pos_EMP-1); %Erase two first reactions: Anab and Decay
aux = 1 - (Kr(pos_EMP:end) == 0); %1s when Kr~=0, and 0s when Kr=0;                                                                          
aux = repmat(aux',size(sto_matrix,1),1); %Create 67 rows of aux'
sto_matrix = sto_matrix.*aux; %Erase those stoichiometric factors for the reactions that not happen???
aux = 1 - (strcmp(pTh.Path, 'Aux')); %1s for non Aux compounds. 0s for Aux compounds
aux = repmat(aux,1, num_reactions); %Create 38 columns of aux.
sto_matrix_eff = ((sto_matrix.*(aux))); %Eliminates from stoMfull the Aux compounds

pos_Met=find(strcmp(Path,'Met')); 
j=1;

for index=1:length(pos_Met)
    logic=sto_matrix_eff(pos_Met(index),:)<0;
    index_nz(index)=sum(logic);
    for i=1:length(logic)
        if logic(i)>0
            matrix(j,i)=1;
            j=j+1;
        else
            matrix(j,i)=0;
        end
    end
end
n_rows=size(matrix,1)-1; %Because last row is added without needing it
matrix=matrix(1:n_rows,:);
RUTA=[zeros(n_rows,(pos_EMP-1)), matrix];

sto_matrix_short=rm.stoMfull;
%sto_matrix_short= (sto_matrix_short>0)-(sto_matrix_short<0);
pos_CO2=strcmp(St.StNames, 'Ci_CO2');
sto_matrix_short(pos_CO2,:) = rm.stoMfull(pos_CO2,:);
sto_matrix_short=sto_matrix_short(1:find(pos_CO2),:);


for i=1:n_rows
    column= RUTA(i,:)>0;
    transp_temp(i,:)=sto_matrix_short(:,column)';
end
no_transport=find(strcmp(Path,'Met'),1,'last'); 
gases=sum(strcmp(St.Phase,'tG'));

%TRANSP=[zeros(n_rows,no_transport),transp_temp(:,no_transport+1:end),zeros(n_rows,gases)];

TRANSP=[zeros(n_rows,no_transport),transp_temp(:,no_transport+1:end)];
pos_Pyr=strmatch('Pyr',rm.rmNames);
num_Pyr=find(RUTA(:,pos_Pyr(end))>0);
TRANSP(1:num_Pyr,:)=TRANSP(1:num_Pyr,:)*2;
TRANSP(num_Pyr+1:end,:)=TRANSP(num_Pyr+1:end,:)*0.5/0.75;


Opt.RUTA=RUTA;
Opt.TRANSP=TRANSP;

%% Inicialization of some variables
Opt.z_r = zeros(rm.num_r,1);
Opt.z_rt = zeros(rmTr.num_tr,1);
Opt.f_r=zeros(rm.num_r,1);

z=zeros(size(RUTA,1),1);
pos_nz=index_nz;
for i=1:length(index_nz)-1
    index_nz(i+1)=index_nz(i+1)+index_nz(i);
end
z(1:index_nz(1))=1/pos_nz(1);

for j=1:length(index_nz)-1
    z(index_nz(j)+1:index_nz(j+1))=1/pos_nz(j+1);
end

%z(index_nz)=1;
Opt.z=z;
Opt.index_nz=index_nz;
Opt.pos_nz=pos_nz;


% Opt.Monod_Glu=(St.Ce_Glu/((Monod(1)+St.Ce_Glu)));
% Opt.Monod_Arg=(St.Ce_Arg/((Monod(2)+St.Ce_Arg)));
% Opt.Monod_Ala=(St.Ce_Ala/((Monod(3)+St.Ce_Ala)));
% Opt.Monod_Asp=(St.Ce_Asp/((Monod(4)+St.Ce_Asp)));
% Opt.Monod_Lys=(St.Ce_Lys/((Monod(5)+St.Ce_Lys)));
% Opt.Monod_Glut=(St.Ce_Glut/((Monod(6)+St.Ce_Glut)));
% Opt.Monod_Ser=(St.Ce_Ser/((Monod(7)+St.Ce_Ser)));
% Opt.Monod_Thr=(St.Ce_Thr/((Monod(8)+St.Ce_Thr)));
% Opt.Monod_Cys=(St.Ce_Cys/((Monod(9)+St.Ce_Cys)));
% Opt.Monod_Gly=(St.Ce_Gly/((Monod(10)+St.Ce_Gly)));
% Opt.Monod_Prol=(St.Ce_Prol/((Monod(11)+St.Ce_Prol)));
% Opt.Monod_Vali=(St.Ce_Vali/((Monod(12)+St.Ce_Vali)));
% Opt.Monod_IsoL=(St.Ce_IsoL/((Monod(13)+St.Ce_IsoL)));
% Opt.Monod_Leu=(St.Ce_Leu/((Monod(14)+St.Ce_Leu)));
% Opt.Monod_Meth=(St.Ce_Meth/((Monod(15)+St.Ce_Meth)));
% Opt.Monod_GluM=(St.Ce_GluM/((Monod(16)+St.Ce_GluM)));
% Opt.Monod_AspG=(St.Ce_AspG/((Monod(17)+St.Ce_AspG)));
% Opt.Monod_Hist=(St.Ce_Hist/((Monod(18)+St.Ce_Hist)));

Opt.pos_Glu=strmatch('Pyr',rm.rmNames);
Opt.pos_Arg=strmatch('Arg',rm.rmNames);
Opt.pos_Ala=strmatch('Ala',rm.rmNames);
Opt.pos_Asp=strmatch('Asp ',rm.rmNames);
Opt.pos_Lys=strmatch('Lys',rm.rmNames);
Opt.pos_Glut=strmatch('Glut',rm.rmNames);
Opt.pos_Ser=strmatch('Ser',rm.rmNames);
Opt.pos_Thr=strmatch('Thr',rm.rmNames);
Opt.pos_Cys=strmatch('Cys',rm.rmNames);
Opt.pos_Gly=strmatch('Gly',rm.rmNames);
Opt.pos_Prol=strmatch('Prol',rm.rmNames);
Opt.pos_Vali=strmatch('Vali',rm.rmNames);
Opt.pos_IsoL=strmatch('IsoL',rm.rmNames);
Opt.pos_Leu=strmatch('Leu',rm.rmNames);
Opt.pos_Meth=strmatch('Meth',rm.rmNames);
Opt.pos_GluM=strmatch('GluM',rm.rmNames);
Opt.pos_AspG=strmatch('AspG',rm.rmNames);
Opt.pos_Hist=strmatch('Hist',rm.rmNames);

pos1 = find(strcmp(St.StNames,'Ce_Glu'));
posend = find(strcmp(St.StNames,'Ce_Hist'));
positions = [pos1 pos1+2:posend];
Opt.Monod_term = St.StV(positions)./(St.StV(positions)+Monod_sat);

Opt.Monod = ones(rm.num_r,1);
Opt.Monod(1:pos_EMP) = 1;

for i=1:St.nSubs
    name=horzcat('Opt.pos_',char(St.SubsNames(i)));
    eval(horzcat('Opt.Monod(',name,')=Opt.Monod_term(i);'));
end
   

R(idR).flagOpt = 0;

%% Output to the global variable R
R(idR).pOp = pOp;         R(idR).pTh = pTh;
R(idR).pKt = pKt;         R(idR).mFd = mFd;      R(idR).Inf_ = Inf_;
R(idR).fullMatrix = fullMatrix;
R(idR).rmTr = rmTr;       R(idR).rm = rm;
R(idR).St = St;
fprintf('DONE!\n')