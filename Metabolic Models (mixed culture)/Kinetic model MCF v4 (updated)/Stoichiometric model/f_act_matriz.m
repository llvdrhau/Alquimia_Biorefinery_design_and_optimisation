function R = f_act_matriz(R,RET)

pos_H2O = strcmp(R.St.StNamesFull,'Ci_H2O');
pos_Pi = find(strcmp(R.St.StNamesFull,'Ci_Pi'));
pos_Pi2 = find(strcmp(R.St.StNames,'Ci_Pi'));

pos_Val = strmatch('Val',R.rm.rmNames);
pos_Ile = strmatch('IsoL',R.rm.rmNames);
pos_Leu = strmatch('Leu',R.rm.rmNames);

%% Valine

R.rm.stoMfull(pos_H2O,pos_Val(1))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Val(1))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Val(1))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Val(1))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Val(1))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Val(2))=RET-2;
R.rm.stoMfull(pos_Pi,pos_Val(2))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Val(2))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Val(2))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Val(2))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Val(3))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Val(3))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Val(3))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Val(3))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Val(3))=3-RET;

R.rm.stoM(pos_Pi2,pos_Val(1):pos_Val(3))=-RET;
R.rm.stoM(pos_Pi2+1,pos_Val(1):pos_Val(3))=RET;
R.rm.stoM(pos_Pi2+2,pos_Val(1):pos_Val(3))=-RET;

%% Isoleucine
R.rm.stoMfull(pos_H2O,pos_Ile(1))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Ile(1))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Ile(1))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Ile(1))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Ile(1))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Ile(2))=RET-2;
R.rm.stoMfull(pos_Pi,pos_Ile(2))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Ile(2))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Ile(2))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Ile(2))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Ile(3))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Ile(3))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Ile(3))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Ile(3))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Ile(3))=3-RET;

R.rm.stoM(pos_Pi2,pos_Ile(1):pos_Ile(3))=-RET;
R.rm.stoM(pos_Pi2+1,pos_Ile(1):pos_Ile(3))=RET;
R.rm.stoM(pos_Pi2+2,pos_Ile(1):pos_Ile(3))=-RET;

%% Leucine
R.rm.stoMfull(pos_H2O,pos_Leu(1))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Leu(1))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Leu(1))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Leu(1))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Leu(1))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Leu(2))=RET-2;
R.rm.stoMfull(pos_Pi,pos_Leu(2))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Leu(2))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Leu(2))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Leu(2))=2-RET;

R.rm.stoMfull(pos_H2O,pos_Leu(3))=RET-3;
R.rm.stoMfull(pos_Pi,pos_Leu(3))=-RET;
R.rm.stoMfull(pos_Pi+1,pos_Leu(3))=RET;
R.rm.stoMfull(pos_Pi+2,pos_Leu(3))=-RET;
R.rm.stoMfull(pos_Pi+3,pos_Leu(3))=3-RET;

R.rm.stoM(pos_Pi2,pos_Leu(1):pos_Leu(3))=-RET;
R.rm.stoM(pos_Pi2+1,pos_Leu(1):pos_Leu(3))=RET;
R.rm.stoM(pos_Pi2+2,pos_Leu(1):pos_Leu(3))=-RET;
end