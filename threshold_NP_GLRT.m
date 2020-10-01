function umbral=threshold_NP_GLRT(Pfa_NP,J,M)

% Determina el umbral para una probabilidad de falsa alarma objetivo

% Numero de muestras
N = 2e6;

% Valores del estadistico
T0 = zeros(N,1);
for n=1:N
    X = chi2rnd(2*M,J,1);
     Y = X - 2*M;
%     Y = X - 2*M*(1+log(X/(2*M)));
     Y=Y.*(Y>0);
    T0(n) = sum(Y);
end


dibuja=0;



%%%%%%%%%%%%%%% Histograma de T(e|H0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dibuja==1
    
figure(30); histogram(T0,'Normalization','pdf') % histograma de T_H0
figure(30);title('Histograma del estadistico T(e|H0)','FontSize',22)
figure(30);xlabel('T','FontSize',20)
figure(30);hold off

end

%%%%%%%%%%%%%%% Pfa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cdf_T0,x] = ecdf(T0);
Pfa = 1-cdf_T0;

indice=find(Pfa <= Pfa_NP,1,'first');
umbral=x(indice);
% disp(' ')
% disp(['umbral: ',num2str(umbral)])

if dibuja==1
    
figure(50);plot(x,cdf_T0,'b',x,Pfa,'r',umbral,Pfa(indice),'ko')
figure(50);title('Pfa(\gamma)','FontSize',22)
figure(50);xlabel('\gamma','FontSize',20)
figure(50);legend('cdf','Pfa','NP threshold')
figure(50);hold off

end



