function I_W = channel_polarization_huawei_approx(N)
% [ref] 3GPP R1-167209 Polar code design and rate matching

    I_W = zeros(1,N);
    C = 2.^((0:log2(N)-1)*0.25);

    for n= 1:N       
        I_W(n) = de2bi(n-1,log2(N),'right-msb') * C';
    end
   
end