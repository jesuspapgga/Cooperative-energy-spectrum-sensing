 function [G]=instantaneous_snrs(Gc,P,Varr)
%
% Calcula la snr en los canales entre los transmisores y los receptores
% suponiendo que solo el transmisor considerado esta activo
%
% Variables de entrada:
% - Gc: Realizaciones de la ganancia en potencia de los canales (N x K x J)
% - P: potencias de transmision (K x 1)
% - Varr: Varianza de ruido en los receptores (J x 1)
%
% Variables de salida:
% - G: SNR's de los canales entre los Tx's y los Rx's (N x K x J)

N=size(Gc,1);
K=size(Gc,2);
if ndims(Gc)==2
    J=1;
else
    J=size(Gc,3); % no. de Rx's
end

G=zeros(N,K,J);

for ntx=1:K
    for nrx=1:J
        G(1:N,ntx,nrx) = P(ntx)*Gc(1:N,ntx,nrx)/Varr(nrx);
    end
end


