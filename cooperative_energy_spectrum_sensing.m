% Centralized cooperative spectrum sensing with soft fusion of energy
% measurements. Neyman-Pearson (NP) detection criterion.
%
%Three methods:
% 
% - LRT (it assumes that the instantaneous SNRs at the cognitive radios (CRs)
%   are known by the fusion center (FC))
%
% - generalized LRT (GLRT)
%
% - Online Expectation-Maximization (EM) based algorithm that jointly estimates the
%   SNRs and detects the presence of primary signals
%


clear all
close all

mySeed=100;
rng(mySeed);

% The sensing channels are assumed to be independent and Rayleigh 
% distributed with the following average SNRs in dBs
SNRdB = [-3,-3,-3];
snr = 10.^(SNRdB/10);

% The sensing channels are time-varying with the following normalized
% Doppler shift (the same for all sensing channels)
fd_n=0.001;

% Number of CRs
J = length(SNRdB);

% Transmit power of the primary user (PU)
P=1;

% Noise variance at the CRs (WGN)
Varr = P./snr;

% No. of sensing samples for energy measurements
M=32; 

% The PU activity is modeled as a Markov chain with the following
% transition probabilities
t00=1/2;
t11=1/2;

% Number of energy measurements
N = 1e3;

% Forgeting factor of the EM online algorithm
mu=0.1; 

% Prescribed probability of false alarm considering the Neyman-Pearson criterion
Pfa_NP=0.01;

% Detectors selection
opt_LRT = 1;
opt_GLRT = 1;
opt_EMO = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multipath fading realizations with time correlation
Gc=sensing_channels_realizations(N,1,J,fd_n,0);

% Instantaneous SNR at the CRs
[G]=instantaneous_snrs(Gc,P,Varr);
g(1:N,1:J) = G(1:N,1,1:J);

% Sequence of PU states (N x 1)
modelot={'markov',t00,t11};
p0 = (1-t11)/(2-t00-t11); % prior probability of H0 
p1=1-p0; % prior probability of H1
s=pu_states(N,modelot); 

% Energy measurements at the CRs location (N x J)
E = energy_measurements(g,s,M); 


%%%%%%%%%%%  LRT detector %%%%%%%%%%%

if opt_LRT == 1
    
% LRT statistic values
T_LRT = LRT_statistics(E,g);

% LRT decision thresholds
metodo = 'aHBE';
umbral_LRT=zeros(N,1);
for n=1:N
    umbral_LRT(n) = umbral_NP_LRT(M,g(n,:)',Pfa_NP,metodo);
end

% Decisions of the LRT detector
se_LRT = zeros(N,1);
for n=1:N
    if T_LRT(n) > umbral_LRT(n)
        se_LRT(n) = 1;
    end
end
Nfa_LRT = sum(se_LRT(s==0)==1);
Nd_LRT = sum(se_LRT(s==1)==1);
Pfa_LRT = Nfa_LRT/sum(s==0);
Pd_LRT = Nd_LRT/sum(s==1);
Pe_LRT = 1-sum(s==se_LRT)/N;

end % opt_LRT

%%%%%%%%%%%  GLRT detector %%%%%%%%%%%

if opt_GLRT == 1

% One-shot Maximum Likelihood (ML) estimates of the SNRs
ge_GLRT = (E/(2*M)-1);
ge_GLRT = ge_GLRT.*(ge_GLRT > 0);

% GLRT statistic values
T_GLRT = LRT_statistics(E,ge_GLRT);

% GLRT decision thresholds 
umbral_GLRT=threshold_NP_GLRT(Pfa_NP,J,M);

% Decisions of the GLRT detector
se_GLRT = zeros(N,1);
for n=1:N
    if T_GLRT(n) > umbral_GLRT
        se_GLRT(n) = 1;
    end
end
Nfa_GLRT = sum(se_GLRT(s==0)==1);
Nd_GLRT = sum(se_GLRT(s==1)==1);
Pfa_GLRT = Nfa_GLRT/sum(s==0);
Pd_GLRT = Nd_GLRT/sum(s==1);
Pe_GLRT = 1-sum(s==se_GLRT)/N;

% Mean square error of the one-shot ML estimates of the SNRs for the GLRT
mse_GLRT=zeros(J,1);
for j=1:J
    mse_GLRT(j)=sum((g(:,j)-ge_GLRT(:,j)).^2)/N;
end
MSE_GLRT = mean(mse_GLRT);

end % opt_GLRT

%%%%%%%%%%%  Online EM detector %%%%%%%%%%%

if opt_EMO == 1
    
% Online EM estimation of the SNRs at the CRs
[ge_EMO,r] = estimates_EM_online(E,M,p0,mu); 

% Values of the EM statistics
T_EMO=zeros(N,1);
T_EMO(1) = sum(E(1,:).*ge_EMO(1,:)./(1+ge_EMO(1,:)));
for n=2:N
    T_EMO(n) = sum(E(n,:).*ge_EMO(n-1,:)./(1+ge_EMO(n-1,:)));
end

% Decision thresholds for the EM detector
umbral_EMO=zeros(N,1); 
metodo = 'aHBE';
umbral_EMO(1) = threshold_NP_GLRT(Pfa_NP,J,M);
for n=2:N
    umbral_EMO(n) = threshold_NP_LRT(M,ge_EMO(n-1,:)',Pfa_NP,metodo);
end

% Decisiones detector EMO
se_EMO = zeros(N,1);
for n=1:N
    if T_EMO(n) > umbral_EMO(n)
        se_EMO(n) = 1;
    end
end
Nfa_EMO = sum(se_EMO(s==0)==1);
Nd_EMO = sum(se_EMO(s==1)==1);
Pfa_EMO = Nfa_EMO/sum(s==0);
Pd_EMO = Nd_EMO/sum(s==1);
Pe_EMO = 1-sum(s==se_EMO)/N;

% MSE of the online EM estimates
mse_EMO=zeros(J,1);
for j=1:J
    mse_EMO(j)=sum((g(:,j)-ge_EMO(:,j)).^2)/N;
end
MSE_EMO = mean(mse_EMO);

end % opt_EMO


% Results
disp(' ')
disp(['Pfa prescribed (NP): ',num2str(Pfa_NP)])
if opt_LRT == 1
    disp(['Pfa_LRT: ',num2str(Pfa_LRT),'    Pd_LRT: ',num2str(Pd_LRT),'     Pe_LRT: ',num2str(Pe_LRT)])
end
if opt_GLRT == 1
    disp(['Pfa_GLRT: ',num2str(Pfa_GLRT),'    Pd_GLRT: ',num2str(Pd_GLRT),'     Pe_GLRT: ',num2str(Pe_GLRT)])
end
if opt_EMO == 1
    disp(['Pfa_EMO: ',num2str(Pfa_EMO),'    Pd_EMO: ',num2str(Pd_EMO),'     Pe_EMO: ',num2str(Pe_EMO)])
end

disp(' ')
if opt_GLRT == 1
    disp(['MSE of the ML estimates for the GLRT: ',num2str(MSE_GLRT)])
end
if opt_EMO == 1
    disp(['MSE of the EM estimates: ',num2str(MSE_EMO)])
end
disp(' ')



