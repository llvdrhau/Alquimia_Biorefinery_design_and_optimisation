function  status = out(t, ~, flag)
%persistent All_StatesVar
status = 0;
switch flag
    case 'init'
       
    case {[]}
        try
            load('All_StatesVar.mat')
        catch
            All_StatesVar=[];
            
        end
        
        R = evalin('base', 'R');
        Opt = evalin('base', 'Opt');
        St = R.St;
        AlgSt = R.AlgSt;        rt = R.rmTr.rt;
        rm = R.rm;              pOp = R.pOp;
        z = Opt.z;
        dX_dt = R.dX_dt;        Qgas = R.Qgas;
        t = t(end);
        
        r=rm.r;
        
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

        nC=[6 3 4 6 5 3 4 3 2 5 5 6 6 5 5 4 6];
        
        r(pos_Arg)=r(pos_Arg)*6;
        r(pos_Ala)=r(pos_Ala)*3;
        r(pos_Asp)=r(pos_Asp)*4;
        r(pos_Lys)=r(pos_Lys)*6;
        r(pos_Glut)=r(pos_Glut)*5;
        r(pos_Ser)=r(pos_Ser)*3;
        r(pos_Thr)=r(pos_Thr)*4;
        r(pos_Cys)=r(pos_Cys)*3;
        r(pos_Gly)=r(pos_Gly)*2;
        r(pos_Prol)=r(pos_Prol)*5;
        r(pos_Vali)=r(pos_Vali)*5;
        r(pos_IsoL)=r(pos_IsoL)*6;
        r(pos_Leu)=r(pos_Leu)*6;
        r(pos_Meth)=r(pos_Meth)*5;
        r(pos_GluM)=r(pos_GluM)*5;
        r(pos_AspG)=r(pos_GluM)*4;
        r(pos_Hist)=r(pos_GluM)*6;
        
         r(pos_Arg)=r(pos_Arg)*176;
        r(pos_Ala)=r(pos_Ala)*96;
        r(pos_Asp)=r(pos_Asp)*96;
        r(pos_Lys)=r(pos_Lys)*224;
        r(pos_Glut)=r(pos_Glut)*144;
        r(pos_Ser)=r(pos_Ser)*80;
        r(pos_Thr)=r(pos_Thr)*128;
        r(pos_Cys)=r(pos_Cys)*144;
        r(pos_Gly)=r(pos_Gly)*48;
        r(pos_Prol)=r(pos_Prol)*176;
        r(pos_Vali)=r(pos_Vali)*192;
        r(pos_IsoL)=r(pos_IsoL)*240;
        r(pos_Leu)=r(pos_Leu)*240;
        r(pos_Meth)=r(pos_Meth)*240;
        r(pos_GluM)=r(pos_GluM)*144;
        r(pos_AspG)=r(pos_AspG)*96;
        r(pos_Hist)=r(pos_Hist)*160;
        
        pos_Ana_Pro=strcmp(rm.rmNames,'AnabProt');
        pos_BM = strcmp(St.StNames,'X');
        pos_EMP = find(strcmp(rm.rmNames,'EMP'));
        pos_Pyr = find(strcmp(rm.rmNames,'Pyr > Ac(For)'));
        
        Ana = r(pos_Ana_Pro)*rm.stoM(pos_BM,pos_Ana_Pro);                  % Production of BM in C-mol BM/Lx·h (Anabolism rate * stoichiometric coeff for BM)
        Prot_Cat = sum(r(pos_EMP+1:pos_Pyr-1));                            % Consumption of AA in catabolism in C-mol AA/Lx·h (Conversion to C-mol from mol was already done)
        Prot_Ana = r(pos_Ana_Pro)*sum(nC)/length(nC);                      % Consumption of AA in anabolism in C-mol AA/Lx·h (anabolism rate * average carbon content of AA)
        
        Yield_PRO= Ana / (Prot_Cat + Prot_Ana);
                
        All_StatesVar(end+1, :) = [t; St.StV; St.Ci_H; St.Ce_H; AlgSt.DGr; rm.r; rt; pOp.Vx; AlgSt.Pcell; Qgas; dX_dt; z; Yield_PRO]';
        fprintf('\n> Current simulation time is %e hours>>\n\n',t);
        assignin('base', 'All_StatesVar', All_StatesVar)
        save All_StatesVar.mat All_StatesVar
        save R.mat R
        save Opt.mat Opt
    case 'done'
        R = evalin('base', 'R');
        Opt = evalin('base', 'Opt');
        %save All_StatesVar.mat All_StatesVar
        save R.mat R
        save Opt.mat Opt
        
        toc
end