function [FRF,FBB] = AltMinPE(Fopt,NT_RF,NT,NS)
ANG =  -pi + (2*pi)*rand(NT,NT_RF);
FRF = exp(1i*ANG);
L = 1:100;
for i = 1:length(L)
    [U,~,V] = svd(Fopt'*FRF);
    V1 = V(:,1:NS);
    FDD = V1*U';
    ANG = angle(Fopt*FDD');
    FRF = exp(1i*ANG);
    FBB = sqrt(NS)*FDD/norm(FRF*FDD,'fro');
    NORM(i) = norm(Fopt - FRF*FBB,inf)^2;
end
    FBB = sqrt(NS)*FDD/norm(FRF*FDD,'fro');
  %  plot(L,NORM);
end