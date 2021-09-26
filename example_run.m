% This test case tests the cdecoding time complexity

clear all;
addpath('Functions');
tStart = tic;

M = 4000; %Codeword length
Rate = 1/3; % coding rate
K = ceil(M*Rate); % infomation length including CRC

%% code construction
isBitRev = false; % wheather bitreversal permutation is used
methods = ["Huawei Approx","heuristic","gaussian approx"];
method = methods(3) ;
designpara = -2; % This parameter is not needed for method(1); 
pcp = PolarCodeStruct(M,K,isBitRev,method,designpara); 

%% for crc aided scl decoder
CRCTypes =  ["none","CRC11"];
CRCType = CRCTypes(1);
nL = 4; % number of list decoders
CRC = crc_generator(CRCType); % Structure CRC
A = K - CRC.L; % infomation bit length 


%% simulation running paramter
EsNodB = -5:0.5:-2;
EsN0 = 10.^(EsNodB/10);
NoVec = 1./EsN0; % Es is normalized to 1
sigma = sqrt(NoVec/2);
[nBitErrs,nBlkErrs,BER_sim,FER_sim] =deal( zeros(size(EsNodB))); 

%% simulation run
t = 0;
for i = 1: length(EsNodB)% loop each SNR(EsN0)
    
    % print progress
    fprintf('\n Now running EbN0 %.1f dB [%d of %d]\n',EsNodB(i),i, length(EsNodB) );
    printLen = 0;
    tic;
    
    % loop each block
    Nblocks = 0;
    while Nblocks < 10000 && nBlkErrs(i)< 100 % stop criterion 
        
        %generate random A-bit message
        msg = randi([0 1],1,A); 
        
        % add CRC
        if CRC.L > 1
            msgwithcrc = [msg,compute_crc(msg,CRC.g)];
        else
            msgwithcrc = msg;
        end

        % encoding
        cword = polar_code_encoder(msgwithcrc,pcp);
        
        % rate matching
        rcw = polar_code_rate_matching(cword,pcp.P);
        
        % BPSK bit to symbol mapping
        s = 1 - 2 * rcw; 
        
        %AWGN channel
        r = s + sigma(i) * randn(1,M); 
        
        % channel llr
        llr = 4 * EsN0(i) *r;
        
        % rate recovery
        llr = polar_code_rate_recovery(llr,pcp.N,pcp.P);
        
        % decding 
        if CRC.L > 1 % SCL decoding
            msg_cap = polar_code_sscl_decoder(llr,pcp,CRC.g,nL);
        else
            msg_cap = polar_code_sc_decoder(llr,pcp);
        end

        %Counting errors
        Nerrs = sum(msg ~= msg_cap);

        if Nerrs > 0
            nBitErrs(i) = nBitErrs(i) + Nerrs;
            nBlkErrs(i) = nBlkErrs(i) + 1;
        end

        Nblocks = Nblocks + 1;
        BER_sim(i) = nBitErrs(i)/K/Nblocks;
        FER_sim(i) = nBlkErrs(i)/Nblocks;

        %print progress
        if mod(Nblocks,10) == 0 
            t = toc; 
            fprintf(repmat('\b',1,printLen));
            printStr = sprintf(' Elasped time is %.1f seconds,# Tx blocks: %d, BER: %.5f, BLER: %.5f',...
                                                                    t,Nblocks,BER_sim(i),FER_sim(i));
            fprintf(printStr);
            printLen  = length(printStr);        
        end
    end


end
    
% display and show simulation results    
fprintf('\n\n');
disp('simlulation results:');
disp('----------------------');
sim.N = pcp.N;
sim.K = pcp.K;
sim.R = Rate;
if pcp.isBitRev
    sim.BitRev = 'yes';
else
    sim.BitRev = 'no';
end
sim.CRCType = CRCType;
if A < pcp.K
    sim.decType = 'sscl';
else
    sim.decType = 'ssc';
end
sim.SNRdes = designpara;
sim.constrctMethod = method;
sim.EbNodB = EsNodB;
sim.BER = BER_sim;
sim.nBitErrs = nBitErrs;
sim.FER = FER_sim;
sim.nBlkErrs = nBlkErrs;
disp(sim);
fileName = sprintf('sim_result_%s_desSNR_%.2fdB_n_15.mat',method,designpara);
save(fileName,'sim');
fprintf('result is saved to file: %s.mat\n',fileName);
fprintf('Total elapsed time is %.2f seconds\n',toc(tStart));


