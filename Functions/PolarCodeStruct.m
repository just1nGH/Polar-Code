% 
classdef PolarCodeStruct < handle
    properties
        M; % the number of coded bits to be transmitted
        K; % the number of information bits and crc parity
        N; % the number of encoded bits as the ouput of polar code encoder
        P; % puncture positions in a vector: 1 by N-M
        U; % frozen bits position in a vector: 1 by N-K
        isBitRev; % indicate if bit reversal permulation is used 
        constructMethod; % the channel polarization methods: 'BP', 'GA','Huawer approx'
        constructDesignPara; % design parameters(not applicable for Huawei approx)
    end
    
    methods
        function obj = PolarCodeStruct(M,K,varargin)
            % M, codeword length
            % K, msg length
            % varargin{1}, if bitreveral permutaion is applied
            % varargin{2},CRC:
            % 'CRC24A','CRC24B','CRC24C','CRC16','CRC11','CRC6'
            
            obj.M = M;
            obj.K = K;
            obj.N = 2^(ceil(log2(M)));
            obj.isBitRev = false;
            obj.constructMethod = 'Huawei Approx';
            obj.constructDesignPara = 0;
            
            if nargin > 5
                    error('too many input parameters!')
            elseif nargin < 2
                    error('you need at lease two inpur arguments!')
            end
            
               
            if nargin >= 3
                obj.isBitRev = varargin{1};
                if nargin >= 4
                    obj.constructMethod = varargin{2};
                    obj.constructDesignPara = 0.01;
                end
                if nargin >= 5
                    obj.constructDesignPara =  varargin{3};
                end
                
            end
                    

            [obj.U,obj.P] = obj.construction();

        end
        
        function [U,P] = construction(obj)
          
           % compute relaible sequence through channel polarization 
           addpath('Functions/channel polarization'); 
           switch lower(obj.constructMethod)
               case 'huawei approx'
                   I_W = channel_polarization_huawei_approx(obj.N);
               case 'heuristic'
                   I_W = channel_polarization_BP(obj.N,obj.constructDesignPara);
               case 'gaussian approx'
                   I_W = channel_polarization_GA(obj.N,obj.constructDesignPara);
               case 'montacalo'
                    load('C:\Simulator\JM0099\01 Simulation Tasks\Task_2_Polar_Code\My Polar Code 2\Functions\channel polarization\weight_MC_N_4096_design_SNR_-3.00.mat','I_W');
               otherwise
                   error('Constrct method is not found!');
           end
           [~,Q] = sort(I_W); Q = Q-1; 
            
            if obj.N == obj.M
               P = [];
               U = Q(1:obj.M-obj.K); 
            else
                
                if obj.isBitRev % with bit reverdal perbutation
                    % Punctured positions
                    %P = bitrevorder(0:obj.N-1); P = P(obj.M+1:obj.N);
                    P = bi2de(de2bi(1:obj.N-obj.M,log2(obj.N),'left-msb'),'right-msb') ; P =P.';
                    % Generate a frozen-set F by selecting the first M-k entries from Q
                    % whose values are smaller than N and do not exist in P
                    Q_hat = setdiff(Q,P,'stable');
                    F = Q_hat(1:obj.M-obj.K);

%                   % Frozen positions
                    U = union(P,F);
                                           
                else% no bit reversal permutation
                    
                    P = Q(1:(obj.N-obj.M));
                    U = Q(1:(obj.N-obj.K));      
                     
                end
            end
            
            P = sort(P);
            U = sort(U);

        end
       
    end
  
end


