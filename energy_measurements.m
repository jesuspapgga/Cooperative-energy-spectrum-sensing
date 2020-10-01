function E=energy_measurements(g,s,M)
%
% Calcula las estimas de energia en cada receptor. 
% Modelo de [Ma2008a]
%
% Variables de entrada:
% - g: SNRs en los receptores (N x J)
% - s: Estados del transmisor (N x 1)
% - M: numero de muestras de los detectores de energia
%
% Variables de salida:
% - E: vectores de energia (N x J)


[N,J] = size(g);

E=zeros(N,J);

for n=1:N
    E(n,1:J) = (s(n)*g(n,1:J)'+1).*chi2rnd(2*M,J,1);
end
