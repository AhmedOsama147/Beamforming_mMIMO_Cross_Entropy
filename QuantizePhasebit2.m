function qFRF = QuantizePhasebit2(FRF,amp)
    ang = angle(FRF);
    ang_quan = zeros(size(ang));
    codebook = pi*[1/4  3/4 -3/4 -1/4];
     qFRF = zeros(size(ang));
    for i = 1:size(ang,1)
        for j = 1:size(ang,2)
            if real(FRF(i,j)) ~= 0
                [~,ind] = min(abs(ang(i,j)-codebook)) ;
                ang_quan(i,j) = codebook(ind);
                 qFRF(i,j) = 1/sqrt(amp)*exp(1i*ang_quan(i,j));
            end
        end
    end
end