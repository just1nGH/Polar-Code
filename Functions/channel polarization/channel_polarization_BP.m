function [I_W] = channel_polarization_BP(N,erasureProb)
% [Ref] Arikan, Erdal. "Channel polarization: A method for constructing 
% capacity-achieving codes for symmetric binary-input memoryless channels." 
% IEEE Transactions on information Theory 55.7 (2009): 3051-3073.

    I_W = zeros(1,2*N-1);

    I_W(1) = 1- erasureProb;

    for i = 2:2:2*N-1
        j = floor(i/2);
        I_W(i) = I_W(j)^2;  
        I_W(i+1) = 2*I_W(j) - I_W(j)^2;
    end
    I_W = I_W(N:end);
    
    
end
