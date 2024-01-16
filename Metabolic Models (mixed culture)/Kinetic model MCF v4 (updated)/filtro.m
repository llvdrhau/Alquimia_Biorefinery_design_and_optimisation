close all
% VFA = iVal;
% for k=1:9
% matriz  = VFA;
% x = k/10;
% y = 1-x;
% for i=2:18;
%     for j=2:20;
%         matriz(i,j) = x*VFA(i,j) + y/4*(VFA(i-1,j) + VFA(i+1,j) + VFA(i,j-1) + VFA(i,j+1));
%     end
% end
% eval(horzcat('matrix_',num2str(k),'=matriz;'))
% f_contorno(matriz)
% end

nomes = {'Ac' 'Pro' 'iBut' 'nBut' 'iVal' 'nVal' 'iCap' 'Et'};

for m=1:length(nomes)
    eval(horzcat('load(''',char(nomes(m)),'_2.mat'')'))
    eval(horzcat('VFA=',char(nomes(m)),';'))
    matriz = VFA;
    k = 5;
    x = k/10;
    y = 1-x;
    for i=2:18;
        for j=2:19;
            matriz(i,j) = x*VFA(i,j) + y/4*(VFA(i-1,j) + VFA(i+1,j) + VFA(i,j-1) + VFA(i,j+1));
        end
    end
    eval(horzcat('matrix_',char(nomes(m)),'=matriz;'))
    f_contorno(matriz)
%     eval(horzcat('save matrix_',char(nomes(m)),'.mat matrix_',char(nomes(m))))
end