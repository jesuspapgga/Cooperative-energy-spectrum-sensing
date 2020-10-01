function [ge,r] = estimates_EM_online(E,M,p0,mu)

% Estima de g con el algoritmo EM online
%
% Variables de entrada:
% - E: matriz de observaciones  (N x J)
% - M: numero de muestras para la estima de la energia
% - p0: probabilidad a priori de PU inactivo
% - mu: factor de olvido
%
% Varibles de salida:
% - ge: estimas de g(n)  (N x J)
% - r: estimas de la probabilidad a posteriori de s(n)=1  (N x 1)


% Numero de observaciones y de sensores
[N,J]=size(E);

p1=1-p0; % probabilidad a priori de PU activo


% Inicializo
gen = E(1,:)'/(2*M)-1;
gen = gen.*(gen>0);
a=zeros(J,1);
b=0;
ge = zeros(N,J);
r = zeros(N,1);

for n=1:N
    
    % E-step
    f01 = exp(log(p0/p1)+sum(M*log(1+gen)-(1/2)*E(n,:)'.*gen./(1+gen)));
    gam = 1/(1+f01);
    r(n) = gam;

    % M-step
    b = (1-mu)*b + mu*2*M*gam;
    a = (1-mu)*a + mu*gam*E(n,:)';
    
    gen = a/b -1;
    gen = gen .* (gen>0);
    
    ge(n,1:J)=gen;
    
end % ciclo a los vectores de observacion
    