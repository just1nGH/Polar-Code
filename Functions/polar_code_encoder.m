function enc = polar_code_encoder(msg,PCP)
%---------------------------------------------------------------------
% Input
%       msg: 1-by-K vector
%       N: length of the codeword 
%       F: Frozen bit postions
%       isRev: is bit reversal to be used ? 0: unused, 1 used. 
%      
% Output
%       codeword: 1-by-M vector
% Juquan Mao 2021 Sep @ juquan.justin.mao@gmail.com
%----------------------------------------------------------------------
    
    % bit-reversal permutation 
    % 3gpp doesn't need this inverse, others do need.
    isBitRev = PCP.isBitRev;
    
    N = PCP.N;
    F = PCP.U;
    
    % put message in the message postion of the codeword
    u = zeros(1,N);
    u(setdiff(0:N-1,F,'stable')+1) = msg; 
 
    enc =  polar_code_encode(u, isBitRev);
   
end


function x = polar_code_encode(u, isBitRev)
    % This function encodes u into x.   x = u*B_N*F_N(isBitRev)   
    % or  x = u*F_N (~isBitRev)
    % where N = 2^n bits, B_N is bit reversal permutaion matrix, F_N = F
    % kroncker N, F = [ 1 0 ; 1 1];
  
    N = size(u,2);
    if (N == 1)
        x = u;
    else
        if isBitRev
            u1_xor_u2 = mod(u(:,1:2:end) + u(:,2:2:end),2);
            u2 = u(2:2:end);
        else
            u1_xor_u2 = mod(u(:,1:N/2) + u(:,N/2+1: N),2);
            u2 = u(:,N/2+1: N);
        end
            
        x = [polar_code_encode(u1_xor_u2,isBitRev) polar_code_encode(u2,isBitRev)];
    end

end
 




 

