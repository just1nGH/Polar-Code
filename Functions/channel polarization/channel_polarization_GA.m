function meanLLR = channel_polarization_GA(N,desSNR)
% N: number of channels must be in power of 2
% desSNR: in dB
% meanLLR: the mean LLR of each polarized channel
% Adapted based on the Ref: H. Ochiai, P. Mitran and H. Vincent Poor, "Capacity-Approaching Polar 
% Codes With Long Codewords and Successive Cancellation Decoding Based on 
% Improved Gaussian Approximation," in IEEE Transactions on Communications, 
% vol. 69, no. 1, pp. 31-43, Jan. 2021
  
    desSNR = 10^(desSNR/10);
    
    gamma = zeros(1,2*N-1);
    gamma(1) = 4*desSNR;
    
    for i = 2:2:2*N-1
        j = floor(i/2);
        gamma(i) = compute_llr_worse_channel(gamma(j));
        gamma(i+1) = compute_llr_better_channel(gamma(j));
    end
    
    meanLLR = gamma(N:end);

end



function gamma = compute_llr_worse_channel(gamma)
 
    u = gamma;

    if u <= 0.2
        gamma = 0.5*u^2 - 0.5*u^3 + (2/3) * u^4;
    else
        z = xi_func(u);
        z = z + log( 2 - exp(z));
        gamma = inv_xi_func(z);        
    end
        
end

function gamma = compute_llr_better_channel(gamma)
    gamma = 2*gamma;
end

function z = xi_func(u)

    a0 = -0.002706; a1 = -0.476711; a2 = 0.0512;
    a = -0.4527; b = 0.0218; c = 0.86;
    Gamma0  = 0.2; Gamma1  = 0.7; Gamma2  = 10; 
    kappa0 = 8.554;

    if u <= Gamma0
        z= -0.5*u + 0.125*u^2 -0.125 * u^3;
    elseif u <= Gamma1
        z = a0 + a1*u + a2 * u^2;  
    elseif u < Gamma2
        z = a* u^c+ b;
    else
        z = -0.25*u + 0.5* log(pi) - 0.5*log(u) +...
            log(1- pi.^2/(4*u) + kappa0/(u^2));
    end

end

function x = inv_xi_func(z)

    Gamma0  = 0.2; Gamma1  = 0.7; Gamma2  = 10; 
    a0 = -0.002706; a1 = -0.476711; a2 = 0.0512;
    a = -0.4527; b = 0.0218; c = 0.86;  
    z0 = xi_func(Gamma0);
    z1 = xi_func(Gamma1);
    z2 = xi_func(Gamma2);
    
    if z >= z0
        x = -2*z + z^2 + z^3;
    elseif z >= z1
        x = (-a1 -  sqrt(a1^2 - 4 * a2*(a0 - z)))/(2 * a2);
    elseif z > z2
        x = ((z- b)/a)^(1/c);
    else
        x = inverse_xi_helper(z);
    end
    
end


function x = inverse_xi_helper(z)

% find reverse function using the bisection search method
    f = @(u)  -0.25*u + 0.5* log(pi) - 0.5*log(u) + log(1- pi/(4*u) + 8.554/(u^2));

    x_l = 10;
    x_r = 100;

    while f(x_r) > z
        x_r = 2*x_r;
    end

    a = f(x_l); b = f(x_r);
    while (a-b) > 1e-1

        x_mid = (x_l + x_r)/2;

        if f(x_mid) < z
            x_r = x_mid;
        else
            x_l = x_mid;
        end

         a = f(x_l); b = f(x_r);
    end

    x = x_mid;

end
