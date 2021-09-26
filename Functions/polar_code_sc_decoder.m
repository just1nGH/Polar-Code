function msgCap = polar_code_sc_decoder(llr,pcp) 
    
    % retreive parameters from pcp struct 
    F = pcp.U;
    N = pcp.N;
    isBitRev = pcp.isBitRev;

    % make F into a bitmap version
    F_bitmap = zeros(1,N); F_bitmap(F+1) = 1;
   
    % decoding
    uCap = node_operations(llr,F_bitmap,isBitRev);
    
    % retrieve the decoded message bits from uCap
    msgCap = uCap(setdiff(0:N-1,F,'stable')+1); 
    
end

function [uCap,beta] = node_operations(alpha,F,isBitRev)
%------------------------------------------------------------------------
% This function implement polor decoding in a recursive manner and utilize
% a binary tree to simplify the decoding process. 'ssc' - simplified
% successive cancellation. 'recur'- recursive
% input
%       alpha: incoing llr of a node, 1-by-N vector 
%       F: Frozen bit position
%       isBitRev: indicate wheather the bitreversal permutaion is used in encoder      
% output
%       beta: hard decision of alpha 1-by-N vector
% Juquan Mao 2021 Sep @ juquan.justin.mao@gmail.com
%------------------------------------------------------------------------

    N = length(alpha);
    
    % inline function f(chck node function) and g (bit node function )
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); %minsum
    % f  = 2* atanh(tanh(0.5*a) .*tanh(0.5*b)); for single parity check code (2,1)
    g = @(a,b,c) b+(1-2*c).*a; %g function for repetition code (2,1)

    if N == 1 % leaf node
        
        if F == 1 % Frozen bit decison is always 0
            beta = 0;
        else      % non-frozen bit decison 0 if llr < 0 , otherwise 1
            beta = (alpha < 0);
        end
        
        uCap = beta; % decoded bit
            
    else % non-leaf node
        
        if isBitRev
             a = alpha(1:2:end);   % 1st half of incoming belief
             b = alpha(2:2:end); % 2nd half of incoming belief 
        else
             a = alpha(1:N/2);   % 1st half of incoming belief
             b = alpha(N/2+1:N); % 2nd half of incoming belief 
        end

         % left handling then move to L child
         [uCapLeft, betaLeft] = node_operations(f(a,b),F(1:N/2),isBitRev); 
         
         % right handling then move to R child  
         c = betaLeft;
         [uCapRight, betaRight] = node_operations(g(a,b,c),F(N/2+1:N),isBitRev);
  
         % U handling
         % compute the hard decision of the node
         beta = [mod(betaLeft + betaRight,2),betaRight];
         if isBitRev
             beta([1:2:end 2:2:end]) = beta;
         end
         
         % concate the decoded bits from its left and right child
         uCap = [uCapLeft uCapRight];

    end
    
    
end
