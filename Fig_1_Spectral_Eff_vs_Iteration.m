% Machine Learning-Based Hybrid Precoding With Low-Resolution Analog Phase Shifters
% IEEE Communication Letters - Jan 2021 - Fig 2
clc; close all; clear;
Nray = 10;
Ncl = 5;
maxCap = 15;
CapStep = 2;
minCap = 0;
minSNR = 0;
maxSNR = 0;
steps = 5;
SNR = minSNR:steps:maxSNR;
EN = 10.^(SNR/10);
varn = 10.^(-SNR/10);
lenSNR = length(SNR);
  
NT_RF  = [2 4 8]; %RF chain at TX
NR_RF  = 4; %RF chain at TX
NT = 64; %Antenna elements at TX
NR = 16; %Antenna elements at RX
NS = 2; %Data Streams
LStyle = ['-.'; '--'];
num = 500;
DEG = 10;
sig = deg2rad(10); % this is the s.d. of the Laplace distribution for the scattering clusters
B = [1 2];
T = [1 5:5:25];
lenT = length(T);
Capacity_CE2 = zeros(num,1);
cap_CE2 = zeros(lenSNR,lenT);

for k = 1:length(NT_RF)
    for i = 1:num
    [H_n,~,~] = mmWaveCh(NT,NR,Ncl,Nray,DEG);
    [U,S,V] = svd(H_n);
   
     Fopt = V(:,1:NS);
     for tT = 1:lenT
         K = 100;
         [FBB_CE_opt,FRF_CE_opt] = cross_entropy_sub2(Fopt,K,T(tT),B(2),NS,NT_RF(k),NT);
         Capacity_CE2(i,tT) = log2(det(eye(NR)+(EN/NS)*H_n*FRF_CE_opt*FBB_CE_opt*FBB_CE_opt'*FRF_CE_opt'*H_n'));   
     end
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for tT = 1:lenT
        cap_CE2(k,tT) = mean(abs(Capacity_CE2(:,tT)));
     end
end

plot(cap_CE2(1,:),'k','LineWidth',2);
hold on;
plot(cap_CE2(2,:),'r','LineWidth',2);
hold on;
plot(cap_CE2(3,:),'b','LineWidth',2);
hold on;
grid on 
legend('$N_t^{RF} = 8$, $N_t^{RF} = 4$,$N_t^{RF} = 2$','Interpreter','latex','FontSize',18)
xlabel('$T$ iteration','FontSize',18,'Interpreter','latex')
xlabel('Spectral Efficiency $bps/Hz$','FontSize',18,'Interpreter','latex')