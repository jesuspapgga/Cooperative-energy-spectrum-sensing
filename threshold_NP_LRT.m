function [umbral]=threshold_NP_LRT(M,g,Pfa_NP,metodo)


if all(g==0)
    
    umbral = inf;
    
else
    
w = g./(1+g);
%w = w(w~=0); % elimina las componentes con peso nulo

J = length(w);

% Aproximacion Gausiana
if strcmp(metodo,'aG')
    umbral = qfuncinv(Pfa_NP) * sqrt(4*M*sum(w.^2)) + 2*M*sum(w);
end

% Aproximacion de Satterthwaite-Welch
if strcmp(metodo,'aSW')
    
    K = 2*M;
    N = K*J; % no. de v.a's chi-square con 1 grado de libertad 
    d=zeros(N,1);
    for j=1:J
        d(K*(j-1)+1:K*j)=w(j);
    end
    
    k=zeros(3,1);
    k(1) = sum(d);
    k(2) = 2*sum(d.^2);
    
    ksw = k(1)^2/k(2);
    tetasw = k(2)/k(1);
    
    umbral = icdf('Gamma',1-Pfa_NP,ksw,tetasw);

end

% Aproximacion de Hall-Buckley-Eagleson
if strcmp(metodo,'aHBE')
    
    K = 2*M;
    N = K*J; % no. de v.a's chi-square con 1 grado de libertad 
    d=zeros(N,1);
    for j=1:J
        d(K*(j-1)+1:K*j)=w(j);
    end
    
    k=zeros(3,1);
    k(1) = sum(d);
    k(2) = 2*sum(d.^2);
    k(3) = 8*sum(d.^3);
    
    eta = 8*k(2)^3/k(3)^2;
    khbe = eta/2;
    tetahbe = 2;
   
    umbral = k(1) + sqrt(k(2)/(2*eta))*(icdf('Gamma',1-Pfa_NP,khbe,tetahbe) - eta);

end

%umbral = umbral - 2*M*sum(log(1+g));

end