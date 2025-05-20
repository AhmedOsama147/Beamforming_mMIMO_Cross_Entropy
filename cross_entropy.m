function [FBB_CE_opt,FRF_CE_opt] = cross_entropy(Fopt,K,T,B,NS,NRF,NT)
    fmin = zeros(T,1);
    Kelite = 0.2*K;
    alpha = 1;
    REF = 1/sqrt(NT)*exp(1j*pi*[-3 -1 1 3]/4);
    %REF = 1/sqrt(NT)*exp(1j*pi*[-1 1]/2);
    wk = zeros(Kelite,1);
    frf = randi([0 2^B-1],NT,NRF,K);
    Pro = (1/2^B)*ones(NT,NRF,2^B);
    theta = mapping2bitresol(frf);
%     FRF_K = 1/sqrt(NT)*exp(1j*pi*theta/2);
    FRF_K = 1/sqrt(NT)*exp(1j*pi*theta/4);
    FRF_CE = zeros(NT,NRF,Kelite);
    FRF_CE_opt =zeros(NT,NRF);
    FBB_CE = zeros(NRF,NS,K);
    FBB_CE_opt = zeros(NRF,NS);
    Fnorm = zeros(K,1);
    for t = 1:T
        if t ~= 1
            for i = 1:NT
                for j = 1:NRF
                    P = reshape(Pro(i,j,:),1,2^B);
                    prob = normalize(P, 'scale', sum(P));
                    cprob = cumsum([0 prob]);
                    x = rand(1,K);
                    frf = discretize(x, cprob)-1;
                    theta = mapping2bitresol(frf);
                    %theta = mapping1bitresol(frf);
                    FRF_K(i,j,:) = 1/sqrt(NT)*exp(1j*pi*theta/4);
                    %FRF_K(i,j,:) = 1/sqrt(NT)*exp(1j*pi*theta/2);
                end
            end
        end
        for k = 1:K
            FBB_CE(:,:,k) = inv(FRF_K(:,:,k)'*FRF_K(:,:,k))*FRF_K(:,:,k)'*Fopt;
            Fnorm(k) = norm(Fopt-FRF_K(:,:,k)*FBB_CE(:,:,k),'fro')^2;
        end
        [~,ind] = mink(Fnorm,Kelite);
        fmin(t) = min(Fnorm);
        FKeNorm = Fnorm(ind);
        SumKelite = sum(Fnorm(ind));
        for ke = 1:Kelite
            wk(ke) = (1/Kelite)*(SumKelite/FKeNorm(ke))^(1);
            FRF_CE(:,:,ke) = FRF_K(:,:,ind(ke));
        end
        
        for i = 1:NT
            for j = 1:NRF
                hh = 0;
                for s = 1:2^B
                    X = 0;
                    SWk = sum(wk) + (1-alpha)*Pro(i,j,s);
                    for m = 1:Kelite
                        if round(angle(FRF_CE(i,j,m)),4) == round(angle(REF(s)),4)
                            X = X + wk(m);
                            hh = hh + 1;
                        end
                        Pro(i,j,s) = alpha * X /SWk;
                    end
                end
            end
        end
    end
 %   plot(fmin);
    for i = 1:NT
        for j = 1:NRF
            [~,index] = max(Pro(i,j,:));
            theta = mapping2bitresol(index-1);
            FRF_CE_opt(i,j) = 1/sqrt(NT)*exp(1j*pi*theta/4);
        end
    end
    FBB_CE_opt = inv(FRF_CE_opt'*FRF_CE_opt)*FRF_CE_opt'*Fopt;
  %  max
    FBB_CE_opt = sqrt(NS)*FBB_CE_opt/(norm(FRF_CE_opt*FBB_CE_opt,'fro'));
end