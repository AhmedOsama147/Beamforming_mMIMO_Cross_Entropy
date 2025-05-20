function [H_n,At,Ar] = mmWaveCh(NT,NR,Ncl,Nray,DEG)
    gamma = sqrt((NT*NR)/(Ncl*Nray));
    alpha = sqrt(1/(2)) *  (randn(Ncl,Nray) + 1i * randn(Ncl,Nray)) ;
    phiT = (-pi) + ((2*pi)*rand(Ncl,1));
    phiR = (-pi) + ((2*pi)*rand(Ncl,1));
    thetaT = (-pi) + ((2*pi)*rand(Ncl,1));
    thetaR = (-pi) + ((2*pi)*rand(Ncl,1));
    sig = deg2rad(DEG); % this is the s.d. of the Laplace distribution for the scattering clusters
    At = [];
    Ar = [];
    H = zeros(NR,NT);
    for j=1:Ncl
         % The final azimuth angles from each cluster i    
         phi_AOD = randLaplacian(Nray, 1, phiT(j), sig^2);
         phi_AOA = randLaplacian(Nray, 1, phiR(j), sig^2);
         theta_AOD = randLaplacian(Nray, 1, thetaT(j), sig^2);
         theta_AOA = randLaplacian(Nray, 1, thetaR(j), sig^2);
         for l=1:Nray
             at = []; ar = [];
             for m=0:1:(sqrt(NT)-1)
                 at = [at;exp(-1j*pi*(m*ones*sin(phi_AOD(l))*sin(theta_AOD(l))+(0:1:sqrt(NT)-1)*cos(theta_AOD(l))))'];
             end
             for m=0:1:(sqrt(NR)-1)
                 ar = [ar;exp(-1j*pi*(m*ones*sin(phi_AOA(l))*sin(theta_AOA(l))+(0:1:sqrt(NR)-1)*cos(theta_AOA(l))))'];
             end
             Alpha = alpha(j,l);
             at = (1/sqrt(NT))*at;
             ar = (1/sqrt(NR))*ar;
             H = H + Alpha*ar*at'; % channel
             At = [At at];
             Ar = [Ar ar];
         end
     end
     H_n = gamma*H;
end