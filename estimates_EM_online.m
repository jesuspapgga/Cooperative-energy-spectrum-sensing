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

ge=0; r=0;    