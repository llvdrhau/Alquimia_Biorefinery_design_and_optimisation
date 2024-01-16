function my_z(R, t)
%% function my_z
%my_z(Opt, R, t)
% returns z vector with the reaction optimisation selection
% Inputs:   Opt         Additional parameter structure
%           t           time
%           R           Parameter structure
Opt = evalin('base', 'Opt');
RUTA=Opt.RUTA;
TRANSP=Opt.TRANSP;
Monod_sat = R.pKt.Monod_sat;
pos_Glu=Opt.pos_Glu;
pos_Arg=Opt.pos_Arg;
pos_Ala=Opt.pos_Ala;
pos_Asp=Opt.pos_Asp;
pos_Lys=Opt.pos_Lys;
pos_Glut=Opt.pos_Glut;
pos_Ser=Opt.pos_Ser;
pos_Thr=Opt.pos_Thr;
pos_Cys=Opt.pos_Cys;
pos_Gly=Opt.pos_Gly;
pos_Prol=Opt.pos_Prol;
pos_Vali=Opt.pos_Vali;
pos_IsoL=Opt.pos_IsoL;
pos_Leu=Opt.pos_Leu;
pos_Meth=Opt.pos_Meth;
pos_GluM=Opt.pos_GluM;
pos_AspG=Opt.pos_AspG;
pos_Hist=Opt.pos_Hist;

pos1 = find(strcmp(R.St.StNames,'Ce_Glu'));
posend = find(strcmp(R.St.StNames,'Ce_Hist'));
positions = [pos1 pos1+2:posend];
Monod_term = R.St.StV(positions)./(R.St.StV(positions)+Monod_sat); %#ok<NASGU>
pos_EMP = find(strcmp(R.rm.rmNames,'EMP'));

Monod = ones(R.rm.num_r,1);

for i=1:R.St.nSubs
    name=horzcat('pos_',char(R.St.SubsNames(i)));
    eval(horzcat('Monod(',name,')=Monod_term(i);'));
end

control=max(abs(Opt.Monod_term-Monod_term));
% Opt.fallo(size(Opt.fallo,1)+1,1)=t;
% Opt.fallo(size(Opt.fallo,1),2:19)=Monod_term-Opt.Monod_term;
% assignin('base', 'Opt', Opt)
%if t>= Opt.T0
if t >= Opt.T0 || control>R.pOp.control
    Opt.Monod = Monod;
    Opt.Monod_term = Monod_term;
    [z] = f_my_optimization(R, Opt,t, RUTA, TRANSP);                       % Optimisation function
    Opt.z = z;
    %% Deactivation of the reactions of those substrate that are not fed into the system
    %Otherwise they would be selected by the model and i) may induce problems
    %in ode due to very little derivatives value and ii) the model would
    %include in the active transport maximum rate determination.
    
    f_r=RUTA*R.AlgSt.f.f_r;
    pos1 = find(strcmp(R.St.StNames,'Ce_Glu'));
    posend = find(strcmp(R.St.StNames,'Ce_Hist'));
    feeding = R.St.StV([pos1 pos1+2:posend]);
    no_feeding = find(feeding<1e-8);
    index = Opt.index_nz;
    desact = ones(size(z,1),1);
    for i=1:length(no_feeding)
        if no_feeding(i)==1
            desact(1:index(1)) = 0;
        else
            desact(index(no_feeding(i)-1)+1:index(no_feeding(i))) = 0;
        end
    end
    z = z.*desact;
    
    %% Convert z vector into metabolic and transport optimised reaction vectors
    z_r = RUTA'*(z-1);
    z_r = max(-1,z_r);                                                         % z_r cannot be lower than -1
    z_r = min(0,z_r);                                                          % z_r cannot be higher than 0
    z_rt = TRANSP'*(z.*f_r);
    z_rt = z_rt - 1;
    
    %% Set next optimising time
    
    if max(abs(z_r - Opt.z_r)) < 1e-3 && Opt.step < 1
        Opt.step = Opt.step*2; %Double step
        Opt.step=min([Opt.step 1]);
    elseif max(abs(z_r - Opt.z_r)) > 0.01 && Opt.step > 1e-3
        Opt.step = Opt.step/2; %Divide step by 2
    end
    Opt.step = max([Opt.step 5e-4]);
    
    fprintf('Step is: %f h.\n',Opt.step)
    
    %% Updating variables
    Opt.z_r = z_r;            Opt.z_rt = z_rt;
    Opt.T0 = t + Opt.step;
    
    assignin('base', 'Opt', Opt)
end
end
