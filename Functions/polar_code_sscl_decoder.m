function msgCap = polar_code_sscl_decoder(llr,pcp,crcG,nL)

    % retrieve paremeters from pcp construct 
    N = pcp.N;
    F = pcp.U;
    A = pcp.K - length(crcG) + 1;  
    isBitRev = pcp.isBitRev;
    
    FinBitMap = zeros(1,N);FinBitMap(F+1) = 1;% make F into a bit-mask
    PM = inf * ones(nL,1); PM(1) = 0;    % path metric  
    llr = repmat(llr,nL,1);% Initiate input llr for all the list decoders
    
    % the recurseive node function
    uCapList = sscl_node_operations(llr,FinBitMap,PM,isBitRev);

    % retrieve the decoded message list from uCap
    msgCapList = uCapList(:,setdiff(0:N-1,F,'stable')+1);

    iSel = 1; %candidate codeword to output, initially set to best PM
    for i = 1:nL
        if crc_check(msgCapList(i,:),crcG) %check if CRC passes
            iSel = i;
            break
        end
    end
    
    % get the winner of the list of decoded message 
     msgCap = msgCapList(iSel,1:A);

     
end


function [uCap,beta,sortnPrunePattern,PM] = sscl_node_operations(alpha,F,PM,isBitRev)
%------------------------------------------------------------------------
% This function implement polor sscl decoding in a recursive manner and utilize
% a binary tree to simplify the decoding process. 'sscl' - simplified
% successive cancellation list. 'recur'- recursive.
% input
%       alpha: of a node, 1-by-N vector 
%       F: Frozen bitmap
%       PM: input path metric  
%       isBitRev: indicates wheather the bitreversal permutaion is used in encoder 
% output
%       uCap: stores the decoded bits 1-by-N vector
%       beta: hard decision of llr 1-by-N vector
%       sortnPrunePattern: the sorted decoder indices selected with nL smallest PMs  
%       PM: updated path metric  
% Juquan Mao 2021 Sep @ juquan.justin.mao@gmail.com
%------------------------------------------------------------------------

    [nL,N]= size(alpha);
    
    % inline function f(L handling) and g (R handling)
    % minsum llr approximation for SPC: llr( mod(a + b,2))
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));
   
    % llr approximation for repetation code: llr( mod(a + b,2))
    g = @(a,b,c) b+(1-2*c).*a; %g function
    
    %% leaf node handling
     if N == 1
         if F == 1
            beta = zeros(nL,N);
            uCap =  beta;
            % add penalty for opposing llr 
            % frozen bit = 0, correnponding llr should be > 0 , if it's
            % < 0, it get a penalty.
            PM = PM + sum(abs(alpha).*(alpha < 0),2); 
            sortnPrunePattern = [];
            
         else
            %Generate 2*nL canidate path
            candiDec = [alpha < 0;alpha >= 0]; %the first half agrees, the second half against
            candiPM =  [PM; PM + abs(alpha)]; % the second half got penalized.

            % Choose nL paths with least PM 
            [PM,pos] = mink(candiPM,nL);
            
            % selected decoders used by its parent to sort and prune candidate decoders
            sortnPrunePattern = pos - (pos > nL)*nL;
            
            % update decisons based on selected decoders
            beta = candiDec(pos);
            
            %hard decsions of the survived decoders
            uCap =  beta;
         end
         
    
    else % non-leaf noded 
        
         %% left handling then move to L child
         if isBitRev
             a = alpha(:,1:2:N);   % 1st half of incoming belief
             b = alpha(:,2:2:N);   % 2nd half of incoming belief              
         else
             a = alpha(:,1:N/2);   % 1st half of incoming belief
             b = alpha(:,N/2+1:N); % 2nd half of incoming belief 
         end
         
         [uCapLeft,betaLeft,sortnPrunePatternLeft,PM] =sscl_node_operations(f(a,b),F(1:N/2),PM,isBitRev); 
         
         %% right handling then move to R child  
         % rearrange the incoming llr based on sort and prune pattern received from its left child 
         for i =1: size(sortnPrunePatternLeft,2)
             a = a(sortnPrunePatternLeft(:,i),:);   % 1st half of incoming belief
             b = b(sortnPrunePatternLeft(:,i),:);  % 2nd half of incoming belief 
         end
         c = betaLeft;
         [uCapRight,betaRight,sortnPrunePatternRight,PM] = sscl_node_operations(g(a,b,c),F(N/2+1:N),PM,isBitRev);
  
         % rearrange the already decoded bits based on sort and prune
         % pattern
         for i = 1: size(sortnPrunePatternRight,2)
             betaLeft = betaLeft(sortnPrunePatternRight(:,i),:);
             uCapLeft = uCapLeft(sortnPrunePatternRight(:,i),:);    
         end
        
         
         %% U handling
         % unite the sort and prune pattern from the left and right child
         sortnPrunePattern = [sortnPrunePatternLeft,sortnPrunePatternRight]; 
         
         % hard decisions corresponding to input alpha 
         beta = [mod(betaLeft + betaRight,2),betaRight];
         if isBitRev
             beta(:,[1:2:end 2:2:end]) = beta;
         end
         
         % decoded bits        
         uCap = [uCapLeft,uCapRight];

    end
    
end

function passed = crc_check(msg,crcg)
    [~,r1] = gfdeconv(fliplr(msg),fliplr(crcg));
    passed = (r1 == 0);
end