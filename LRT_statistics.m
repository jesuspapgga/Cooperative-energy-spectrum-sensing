function T = LRT_statistics(E,g)

% Detector LRT. Cooperative detection.
%
% Variables de entrada:
% - E: matriz de observaciones  (N x J)
% - g: SNR primaria (N x J)
%
% Varibles de salida:
% - T: valor del estadistico (N x 1)

N=size(E,1);

T=zeros(N,1);

for n=1:N

    T(n) = sum(E(n,:).*g(n,:)./(1+g(n,:)));

end
 
