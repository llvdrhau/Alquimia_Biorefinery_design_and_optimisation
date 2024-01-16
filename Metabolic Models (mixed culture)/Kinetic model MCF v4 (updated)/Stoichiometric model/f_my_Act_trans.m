function [R] = f_my_Act_trans(R, Opt)
%% function f_my_Act_trans
% [R] = f_my_Act_trans(R, Opt)
% returns structure R updating:
% - active transport rates
% - substrate-related transport rates
% Inputs:   R           parameter structure
%           Opt         additional parameter structure
% Outputs:
%           R           parameter structure

%% Loading variables
St = R.St;
pOp = R.pOp;
rmTr = R.rmTr;
DGtr = R.AlgSt.DGtr;
numrt = length(DGtr); 
K_act = R.pKt.KTr.K_actV(1:numrt); 
M_act = R.pKt.KTr.M_actV(1:numrt); 
StTrans = St.StV(1:numrt); 
rt_diff  = rmTr.rt_dif(1:numrt); 
z_rt = Opt.z_rt(1:numrt); 

%% Active transport rate calculation
rt_act = K_act*St.X.*((StTrans)./((M_act)+StTrans));                       % (mol_i/Lr h)
rt_act = rt_act.*(1 + z_rt);                                               % Only those active transport reaction selected by the optimiser

switchFunction = (tanh(R.pOp.Gan*(StTrans-0.01))+1)/2;                           % Step function

reac=(R.rm.stoM(1:numrt,:)*R.rm.r)/(pOp.Vr/pOp.Vx).*sign(K_act);           % Reaction rates of actively transported compounds (mol_i/Lr h)
diffusion= rt_diff.*sign(K_act);                                           % Diffusion transport rates of actively transported compounds (mol_i/Lr h)
A=reac-diffusion;                                                          
A=max(A,0);  
rt_act = switchFunction.*(A) + (1-switchFunction).*rt_act;                 % When step function is active, active transport is equal to the productoin rate minus the diffusion transport

%% Total transport rate calculation
rt  = rt_diff + rt_act;                                                    % (mol_i/Lr h)
rt = [rt;rmTr.rt_dif(numrt+1:end)];                                            % Add gas compunds transport rates 

%% Substrate transport rate calculation
pos_AnaG=strcmp(R.rm.rmNames,'Anab');
pos_AnaP=strcmp(R.rm.rmNames,'AnabProt');

index_glu = strcmp(rmTr.rmTrNames,'GluTr');
pos_EMP = find(strcmp(R.rm.rmNames,'EMP'));
reac = sum(R.rm.r(pos_EMP));                                               % (mol_i/Lx h)                            
ana = -R.rm.stoM(index_glu,pos_EMP)*R.rm.r(pos_AnaG);                      % (mol_i/Lx h)
rt(index_glu) = -(reac + ana)*(pOp.Vx/pOp.Vr);                             % (mol_i/Lr h)

index_Arg= strcmp(rmTr.rmTrNames,'ArgTr');
pos_Arg=Opt.pos_Arg;
reac = sum(R.rm.r(pos_Arg));
ana = -R.rm.stoM(index_Arg,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Arg) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Ala= strcmp(rmTr.rmTrNames,'AlaTr');
pos_Ala=Opt.pos_Ala;
reac = sum(R.rm.r(pos_Ala));
ana = -R.rm.stoM(index_Ala,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Ala) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Asp= strcmp(rmTr.rmTrNames,'AspTr');
pos_Asp=Opt.pos_Asp;
reac = sum(R.rm.r(pos_Asp));
ana = -R.rm.stoM(index_Asp,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Asp) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Lys= strcmp(rmTr.rmTrNames,'LysTr');
pos_Lys=Opt.pos_Lys;
reac = sum(R.rm.r(pos_Lys));
ana = -R.rm.stoM(index_Lys,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Lys) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Glut= strcmp(rmTr.rmTrNames,'GlutTr');
pos_Glut=Opt.pos_Glut;
reac = sum(R.rm.r(pos_Glut));
ana = -R.rm.stoM(index_Glut,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Glut) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Ser= strcmp(rmTr.rmTrNames,'SerTr');
pos_Ser=Opt.pos_Ser;
reac = sum(R.rm.r(pos_Ser));
ana = -R.rm.stoM(index_Ser,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Ser) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Thr= strcmp(rmTr.rmTrNames,'ThrTr');
pos_Thr=Opt.pos_Thr;
reac = sum(R.rm.r(pos_Thr));
ana = -R.rm.stoM(index_Thr,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Thr) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Cys= strcmp(rmTr.rmTrNames,'CysTr');
pos_Cys=Opt.pos_Cys;
reac = sum(R.rm.r(pos_Cys));
ana = -R.rm.stoM(index_Cys,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Cys) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Gly= strcmp(rmTr.rmTrNames,'GlyTr');
pos_Gly=Opt.pos_Gly;
reac = sum(R.rm.r(pos_Gly));
ana = -R.rm.stoM(index_Gly,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Gly) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Prol= strcmp(rmTr.rmTrNames,'ProlTr');
pos_Prol=Opt.pos_Prol;
reac = sum(R.rm.r(pos_Prol));
ana = -R.rm.stoM(index_Prol,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Prol) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Vali= strcmp(rmTr.rmTrNames,'ValiTr');
pos_Vali=Opt.pos_Vali;
reac = sum(R.rm.r(pos_Vali));
ana = -R.rm.stoM(index_Vali,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Vali) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_IsoL= strcmp(rmTr.rmTrNames,'IsoLTr');
pos_IsoL=Opt.pos_IsoL;
reac = sum(R.rm.r(pos_IsoL));
ana = -R.rm.stoM(index_IsoL,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_IsoL) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Leu= strcmp(rmTr.rmTrNames,'LeuTr');
pos_Leu=Opt.pos_Leu;
reac = sum(R.rm.r(pos_Leu));
ana = -R.rm.stoM(index_Leu,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Leu) = -(reac + ana)*(pOp.Vx/pOp.Vr);   

index_Meth= strcmp(rmTr.rmTrNames,'MethTr');
pos_Meth=Opt.pos_Meth;
reac = sum(R.rm.r(pos_Meth));
ana = -R.rm.stoM(index_Meth,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Meth) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_GluM= strcmp(rmTr.rmTrNames,'GluMTr');
pos_GluM=Opt.pos_GluM;
reac = sum(R.rm.r(pos_GluM));
ana = -R.rm.stoM(index_GluM,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_GluM) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_AspG= strcmp(rmTr.rmTrNames,'AspGTr');
pos_AspG=Opt.pos_AspG;
reac = sum(R.rm.r(pos_AspG));
ana = -R.rm.stoM(index_AspG,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_AspG) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

index_Hist= strcmp(rmTr.rmTrNames,'HistTr');
pos_Hist=Opt.pos_Hist;
reac = sum(R.rm.r(pos_Hist));
ana = -R.rm.stoM(index_Hist,pos_AnaP)*R.rm.r(pos_AnaP);
rt(index_Hist) = -(reac + ana)*(pOp.Vx/pOp.Vr);  

aux = [zeros(numrt,1);ones((length(rmTr.rt_dif)-numrt),1)];                    % Select the liquid components
rt = (1-aux).*rt*(pOp.Vr/pOp.Vliq) + aux.*rt;                              % (mol_i/Lliq h) for all reactions

%% Variable update
R.rmTr.rt_act = rt_act;                                                    % (mol_i/Lr·h)
R.rmTr.rt = rt;                                                            % (mol_i/Lliq·h)
