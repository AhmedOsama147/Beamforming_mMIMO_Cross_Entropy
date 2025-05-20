% Machine Learning-Based Hybrid Precoding With Low-Resolution Analog Phase Shifters
% IEEE Communication Letters - Jan 2021 
clc; close all; clear;
Nray = 10;
Ncl = 5;
maxCap = 15;
CapStep = 2;
minCap = 0;
minNTRF = 1;
maxNTRF = 5;
steps = 1;

SNR = 0;
EN = 10.^(SNR/10);
varn = 10.^(-SNR/10);
lenSNR = length(SNR);
  
NT_RF  = [2 4 8 16 32]; %RF chain at TX
lenNTRF = length(NT_RF);
NR_RF  = 4; %RF chain at TX
NT = 64; %Antenna elements at TX
NR = 16; %Antenna elements at RX
NS = 2; %Data Streams
LStyle = ['-.'; '--'];
num = 500;
DEG = 10;
sig = deg2rad(10); % this is the s.d. of the Laplace distribution for the scattering clusters
B = [1 2];
Capacity_fully_digital = zeros(num,1);
cap_full_digital = zeros(lenSNR,1);
Capacity_PE1 = zeros(num,1);
cap_PE1 = zeros(lenSNR,1);
Capacity_PE2 = zeros(num,1);
cap_PE2 = zeros(lenSNR,1);
Capacity_CE2 = zeros(num,1);
cap_CE2 = zeros(lenSNR,1);
Capacity_OMP2 = zeros(num,1);
cap_OMP2 = zeros(lenSNR,1);

for k = 1:lenNTRF
    for i = 1:num
    [H_n,AT,AR] = mmWaveCh(NT,NR,Ncl,Nray,DEG);
    [U,S,V] = svd(H_n);
   
     Fopt = V(:,1:NS);
     Capacity_fully_digital(i) = log2(det(eye(NR)+(EN/NS)*H_n*Fopt*Fopt'*H_n'));
     
     [FRF_OMP,FBB_OMP] = OMP(Fopt,NT_RF(k),NT/NT_RF(k),AT,NS);
     qFRF_OMP = QuantizePhasebit2(FRF_OMP,NT);
     Capacity_OMP2(i) = log2(det(eye(NR)+(EN/NS)*H_n*qFRF_OMP*FBB_OMP*FBB_OMP'*qFRF_OMP'*H_n'));

     [FRF_PE,FBB_PE] = AltMinPE(Fopt,NT_RF(k),NT,NS);
     qFRF_PE = QuantizePhasebit2(FRF_PE,1);
     Capacity_PE2(i) = log2(det(eye(NR)+(EN/NS)*H_n*qFRF_PE*FBB_PE*FBB_PE'*qFRF_PE'*H_n'));   
%      
     K = 150; T = 40;
     [FBB_CE_opt,FRF_CE_opt] = cross_entropy(Fopt,K,T,B(2),NS,NT_RF(k),NT);
     Capacity_CE2(i) = log2(det(eye(NR)+(EN/NS)*H_n*FRF_CE_opt*FBB_CE_opt*FBB_CE_opt'*FRF_CE_opt'*H_n'));   


     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cap_full_digital(k) = mean(abs(Capacity_fully_digital));
cap_OMP2(k) = mean(abs(Capacity_OMP2));
cap_PE2(k) = mean(abs(Capacity_PE2));
cap_CE2(k) = mean(abs(Capacity_CE2));
end


plot(NT_RF,cap_full_digital(:),"Color",'k','LineWidth',2,'Marker','diamond');
hold on;
plot(NT_RF,cap_CE2(:),"Color",'r','LineStyle','--','LineWidth',2,'Marker','square');
hold on;
plot(NT_RF,cap_OMP2(:),"Color",[0.9 0.5 0.5],'LineStyle',':','LineWidth',2,'Marker','o');
hold on;
plot(NT_RF,cap_PE2(:),"Color",'b','LineStyle','--','LineWidth',2,'Marker','.');
hold on;
legend('Fully Digital','Machine Learning', 'PE-AltMin','OMP Fully Connected',...
   'fontsize',18,'Interpreter','latex')
xlabel('$SNR_{dB}','FontSize',14,'Interpreter','latex')
ylabel('Spectral Efficiency (bps/Hz)','FontSize',14,'Interpreter','latex')