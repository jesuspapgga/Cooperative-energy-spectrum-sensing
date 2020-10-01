function Gc=sensing_channels_realizations(N,K,J,fd_n,Kr)
%
% Genera una realizacion de los canales (Rice) entre los Tx's y los 
% Rx's
%
% Variables de entrada:
% - N: no. de realizaciones de los canales
% - J: no. de Rx's
% - K: no. de Tx's
% - fd_n: maximum normalized Doppler frecuency
% - Kr: factor de Rice
%
% Variables de salida:
% - Gc: realizacion de la ganancia en potencia de los canales (N x K x J)

grafica=0;

Gc=zeros(N,K,J);
for nrx=1:J
    for ntx=1:K
        h=rice_model(Kr,N,fd_n); % (N x 1)        
        g=abs(h).^2; % (N x 1)
        Gc(1:N,ntx,nrx)=g;
    end
end

if grafica == 1
contador=0;
for ntx=1:K
    for nrx=1:J
        contador=contador+1;
        figure(40);subplot(J,K,contador);plot(20*log10(Gc(:,ntx,nrx)),'bo-')
    end
end
figure(40);title('Evolucion de los canales con el tiempo')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=rice_model(Kr,N,fd_n)
% It computes N time-correlated channel gain samples from a fading
% Rice channel with maximum normalized Doppler frequency fd_n. 
% It uses the Jake's model of MATLAB. 

canal = comm.RayleighChannel('SampleRate',1,'MaximumDopplerShift',fd_n);
sig = 1i*ones(N,1); 
h_rayleigh = canal(sig); % (N x 1)

h_LOS=(1+1i)/sqrt(2);

h=sqrt(Kr/(1+Kr))*h_LOS+sqrt(1/(1+Kr))*h_rayleigh; % (N x 1)

end
